#!/usr/bin/env nextflow
// Include processes, alphabetic order of process alias
include { CheckQC } from './CustomModules/CheckQC/CheckQC.nf'
include { EditSummaryFileHappy } from './CustomModules/Utils/EditSummaryFileHappy.nf'
include { GATK4_SELECTVARIANTS } from './modules/nf-core/gatk4/selectvariants/main' 
include { HAPPY_HAPPY as HAPPY_HAPPY_pairwise} from './modules/nf-core/happy/happy/main' 
include { HAPPY_HAPPY as HAPPY_HAPPY_tp_giab} from './modules/nf-core/happy/happy/main' 
include { HAPPY_HAPPY } from './modules/nf-core/happy/happy/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main' 
include { VersionLog } from './CustomModules/Utils/VersionLog.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

def analysis_id = params.outdir.split('/')[-1]

workflow {
    def createMetaWithIdName = {file -> [[id: file.name], file]}
    def addTmpId = {meta, file -> [meta.id, meta, file]}
    def createHappyInput = {meta_truth, truth, meta_query, query -> 
        meta = [
            id: meta_query.id + "_" + meta_truth.id,
            query: meta_query.id,
            truth: meta_truth.id
        ]
        return [meta, query, truth, regions_bed, targets_bed]}

    // Empty channel for optional inputs, where meta val is required (often input tuples)
    empty = Channel.of([[id: "null"], []]).first()

    // Reference file channels
    assembly_to_use = params[params.nist_version_to_use].assembly
    ch_fasta = Channel.fromPath("${params.assembly[assembly_to_use].ref_fasta}").map(createMetaWithIdName).first()
    ch_fasta_fai = Channel.fromPath("${params.assembly[assembly_to_use].ref_fai}").map(createMetaWithIdName).first()
   
    // Reference bed files
    regions_bed = "${params[params.nist_version_to_use].high_conf_bed}"
    targets_bed = "${params.exome_target_bed}"

    // GIAB reference file channels
    ch_giab_truth = Channel.fromPath("${params[params.nist_version_to_use].truth_vcf}")
    .map{file -> 
        tokens = file.name.tokenize("_")
        [[id: tokens[0] + "_truth"], file]
    }
    .first()

    // Input vcf file channel
    ch_vcf_files = Channel.fromPath(["${params.vcf_path}/*.vcf", "${params.vcf_path}/*.vcf.gz"])
    .map { vcf ->
        tokens = vcf.name.tokenize(params.delim)
        id_items = params.id_index.collect{idx -> tokens[idx]}
        identifier = (id_items.join("_")? id_items.join("_") : id_items)
        meta = [
            id: identifier,
            vcf: vcf.simpleName,
            single_end:false
        ]
        [meta, vcf]
    }

    // Get all combinations of unordered vcf pairs, without self-self and where a+b == b+a 
    def lst_used = []

    // Create a channel with all vcf files and combine with input vcf files
    ch_vcf_pairwise = ch_vcf_files.combine(ch_vcf_files)
    .branch {meta_truth, truth, meta_query, query ->
        meta = [
            id: meta_query.id + "_" + meta_truth.id,
            query: meta_query.id,
            truth: meta_truth.id
        ]
        // Select valid combination: is without self-self and only unique sets (a+b == b+a)
        valid: (
                meta_query.id != meta_truth.id
                && lst_used.indexOf(meta_query.id + "_" + meta_truth.id) == -1 
                && lst_used.indexOf(meta_truth.id + "_" + meta_query.id) == -1
            )
            lst_used.add(meta_query.id + "_" + meta_truth.id)
            return [meta, query, truth, regions_bed, targets_bed]
    }

    // Create a channel with vcfs against giab.
    ch_vcf_giab = ch_vcf_files.combine(ch_giab_truth)
    .map(createHappyInput)

    HAPPY_HAPPY(ch_vcf_giab, ch_fasta, ch_fasta_fai, empty, empty, empty)


    // Retrieve true-positives from pairwise comparisons.
    HAPPY_HAPPY_pairwise(ch_vcf_pairwise, ch_fasta, ch_fasta_fai, empty, empty, empty)
    ch_pairwise_vcf_index = HAPPY_HAPPY_pairwise.out.vcf.map(addTmpId)
        .join(HAPPY_HAPPY_pairwise.out.tbi.map(addTmpId), by: 0)
        .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, []]} 
    // SelectVariants on VCF + index to select true-positives
    GATK4_SELECTVARIANTS(
        HAPPY_HAPPY_pairwise.out.vcf.map(addTmpId)
        .join(HAPPY_HAPPY_pairwise.out.tbi.map(addTmpId), by: 0)
        .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, []]}
    )
    // Run HAPPY on pairwise true-positives against GIAB truth
    HAPPY_HAPPY_tp_giab(
        GATK4_SELECTVARIANTS.out.vcf.combine(ch_giab_truth).map(createHappyInput), 
        ch_fasta, ch_fasta_fai, empty, empty, empty
    )

    EditSummaryFileHappy(
        Channel.empty().mix(
            HAPPY_HAPPY.out.summary_csv,
            HAPPY_HAPPY_tp_giab.out.summary_csv,
        )
    )
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
            GATK4_SELECTVARIANTS.out.versions,
            HAPPY_HAPPY.out.versions,
            HAPPY_HAPPY.out.summary_csv.map{meta, csv -> [csv]},
            HAPPY_HAPPY_tp_giab.out.versions,
            HAPPY_HAPPY_tp_giab.out.summary_csv.map{meta, csv -> [csv]},
            CheckQC.out.qc_output,
            VersionLog.out.versions,
            Workflow_ExportParams.out
        ).collect(), 
        multiqc_yaml, [], []
    )
}
