#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CheckQC } from './CustomModules/CheckQC/CheckQC.nf'
include { EditSummaryFileHappy } from './CustomModules/Utils/EditSummaryFileHappy.nf'
include { HAPPY_HAPPY } from './modules/nf-core/happy/happy/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main' 

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
        regions_bed = "/hpc/diaggen/data/databases/GIAB/NIST_v3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
        targets_bed = "/hpc/diaggen/software/production/Dx_tracks/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank.bed"
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
    MULTIQC(analysis_id, HAPPY_HAPPY.out.summary_csv)
}
