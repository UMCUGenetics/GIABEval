#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CheckQC } from './CustomModules/CheckQC/CheckQC.nf'
include { EditSummaryFileHappy } from './CustomModules/Utils/EditSummaryFileHappy.nf'
include { HAPPY_HAPPY } from './modules/nf-core/happy/happy/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main' 
include { VersionLog } from './CustomModules/Utils/VersionLog.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

def analysis_id = params.outdir.split('/')[-1]

workflow {
    def create_meta_with_id_name = {file -> [[id: file.name], file]}
    // reference file channels
    ch_fasta = Channel.fromPath("${params.ref_fasta}").map(create_meta_with_id_name).first()
    ch_fasta_fai = Channel.fromPath("${params.ref_fai}").map(create_meta_with_id_name).first()

    // GIAB reference file channels
    ch_giab_truth = Channel.fromPath("${params["NIST_v2_19"].truth_vcf}").map(create_meta_with_id_name).first()

    // input vcf file channel
    ch_vcf_files = Channel.fromPath(["${params.vcf_path}/*.vcf", "${params.vcf_path}/*.vcf.gz"])
    .map { vcf ->
        (sample, date, flowcell, runnr, barcode, projectname, projectnr) = vcf.name.tokenize("_")
        meta = [
            id: sample,
            single_end:false
        ]
	    [meta, vcf]
    }

    // get all combinations of unordered vcf pairs, without self-self and where a+b = b+a 
    def lst_used = []
    ch_comb = ch_vcf_files.concat(ch_giab_truth).combine(ch_vcf_files)
    .branch {mq, q, mt, t ->
        regions_bed = "${params["NIST_v2_19"].high_conf_bed}"
        targets_bed = "${params.exome_target_bed}"
        meta = [
            id: mq.id + "_" + mt.id,
            query: mq.id,
            truth: mt.id
        ]
        // valid combination is without self-self and only unique sets (a+b = b+a)
        valid: mq.id != mt.id && lst_used.indexOf(mq.id + "_" + mt.id) == -1 && lst_used.indexOf(mt.id + "_" + mq.id) == -1
            lst_used.add(mq.id + "_" + mt.id)
            return [meta, q, t, regions_bed, targets_bed]
    }
    empty = Channel.of([[id: "emptychannel"], []]).first()
    HAPPY_HAPPY(ch_comb, ch_fasta, ch_fasta_fai, empty, empty, empty)
    EditSummaryFileHappy(HAPPY_HAPPY.out.summary_csv)
    CheckQC(
        analysis_id,
        Channel.empty().mix(
            EditSummaryFileHappy.out.indel_all_csv,
            EditSummaryFileHappy.out.snp_all_csv,
        ).collect()
    )

    // Create log files: Repository versions and Workflow params
    VersionLog(Channel.of("${workflow.projectDir}/"))
    Workflow_ExportParams()
    
    multiqc_yaml = Channel.fromPath("${params.multiqc_yaml}")
    MULTIQC(
        Channel.empty().mix(
            HAPPY_HAPPY.out.versions,
            HAPPY_HAPPY.out.summary_csv.map{meta, csv -> [csv]},
            CheckQC.out.qc_output,
            VersionLog.out.versions
        ).collect(), 
        multiqc_yaml, [], [])
}
