#!/usr/bin/env nextflow
// Include processes, alphabetic order of process alias
include { BCFTOOLS_ANNOTATE } from './modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_INPUT } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_GIAB } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_REHEADER } from './modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INPUT } from './modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_GIAB} from './modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_PRIMARY} from './modules/nf-core/bcftools/view/main'
include { CheckQC } from './CustomModules/CheckQC/CheckQC.nf'
include { EditSummaryFileHappy } from './CustomModules/Utils/EditSummaryFileHappy.nf'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_NOCALL } from './modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_TP } from './modules/nf-core/gatk4/selectvariants/main'
include { HAPPY_HAPPY as HAPPY_HAPPY_single } from './modules/nf-core/happy/happy/main'
include { HAPPY_HAPPY as HAPPY_HAPPY_pairwise} from './modules/nf-core/happy/happy/main' 
include { HAPPY_HAPPY as HAPPY_HAPPY_tp_giab} from './modules/nf-core/happy/happy/main' 
include { MULTIQC } from './modules/nf-core/multiqc/main' 
include { VersionLog } from './CustomModules/Utils/VersionLog.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

def analysis_id = params.outdir.split('/')[-1]

// Add banner via log.info
log.info """\
    G I A B E V A L   P I P E L I N E
    ===================================
    input: ${params.vcf_path}
    output: ${params.outdir}
    nist version: ${params.nist_version_to_use}
    assembly: ${params.assembly[params[params.nist_version_to_use].assembly]}
    GIAB settings: ${params[params.nist_version_to_use]}
    ===================================
    
    """
.stripIndent(true)

workflow {
    def createMetaWithIdName = {file -> [[id: file.name], file]}
    def addTmpId = {meta, file -> [meta.id, meta, file]}
    def createHappyInput = {meta_query, query, meta_truth, truth ->
        meta = [
            id: meta_query.id + "_" + meta_truth.id,
            query: meta_query.id,
            truth: meta_truth.id
        ]
        return [meta, query, truth, regions_bed, targets_bed]
    }

    // Empty channel for optional inputs, where meta val is required (often input tuples)
    empty = Channel.of([[id: "null"], []]).first()

    // Reference file channels
    assembly_to_use = params[params.nist_version_to_use].assembly
    ch_fasta = Channel.fromPath("${params.assembly[assembly_to_use].ref_fasta}").map(createMetaWithIdName).first()
    ch_fasta_fai = Channel.fromPath("${params.assembly[assembly_to_use].ref_fai}").map(createMetaWithIdName).first()
   
    // Reference bed files
    regions_bed = "${params[params.nist_version_to_use].high_conf_bed}"
    targets_bed = "${params.assembly[params[params.nist_version_to_use].assembly].exome_target_bed}"

    // GIAB reference file channels
    ch_giab_truth = Channel.fromPath("${params[params.nist_version_to_use].truth_vcf}")
    .map{file ->
        tokens = file.name.tokenize("_")
        [[id: tokens[0] + "_truth"], file]
    }
    .first()

    // Input vcf file channel
    ch_vcf_files = Channel.fromPath(["${params.vcf_path}/*.vcf.gz", "${params.vcf_path}/*.vcf"])
    .map { vcf ->
        // Split filename using params.delim and select indices to create unique identifier
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

    //Compress all VCF and index
    BCFTOOLS_VIEW_INPUT(ch_vcf_files.map { meta, vcf -> [ meta, vcf, [] ]}, Channel.empty().toList(), Channel.empty().toList(), Channel.empty().toList())
    BCFTOOLS_VIEW_GIAB(ch_giab_truth.map { meta, vcf -> [ meta, vcf, [] ]}, Channel.empty().toList(), Channel.empty().toList(), Channel.empty().toList())

    // Slice Input VCFs for primary contigs
    BCFTOOLS_VIEW_PRIMARY(BCFTOOLS_VIEW_INPUT.out.vcf
        .join(BCFTOOLS_VIEW_INPUT.out.tbi)
        .map { meta, vcf, tbi -> [ meta, vcf, tbi ]}, Channel.empty().toList(), Channel.empty().toList(), Channel.empty().toList()
    )

    /*
    BCFTOOLS_NORM (normalisation) is required to
        - place an indel at the left-most position (left-align)
        - normalizes split multiallelic sites into biallelics
    */
    BCFTOOLS_NORM_INPUT(BCFTOOLS_VIEW_PRIMARY.out.vcf.join(BCFTOOLS_VIEW_PRIMARY.out.tbi), ch_fasta)
    BCFTOOLS_NORM_GIAB(ch_giab_truth.join(BCFTOOLS_VIEW_GIAB.out.tbi), ch_fasta)

    // Create a channel with vcfs against giab.
    ch_vcf_giab = BCFTOOLS_NORM_INPUT.out.vcf
    .combine(BCFTOOLS_NORM_GIAB.out.vcf)
    .map(createHappyInput)

    // Get all combinations of unordered vcf pairs, without self-self and where a+b == b+a
    def lst_used = []

    // Create a channel with all vcf files and combine with input vcf files
    ch_vcf_pairwise = BCFTOOLS_NORM_INPUT.out.vcf
    .combine(BCFTOOLS_NORM_INPUT.out.vcf)
    .branch {meta_truth, truth, meta_query, query ->
        meta = [
            id: "pairwise_" + meta_query.id + "_" + meta_truth.id,
            query: meta_query.id,
            truth: meta_truth.id
        ]
        /*
        Select valid combination:
         - without self-self
         - only for sorted pair to ensure correct orientation based on sampleID
           irrespective of input order in channel(s)
        */
        valid: (
            meta_query.id != meta_truth.id
            && [meta_query.id, meta_truth.id] == [meta_query.id, meta_truth.id].sort()
        )
        lst_used.add("pairwise_" + meta_query.id + "_" + meta_truth.id)
        return [meta, query, truth, regions_bed, targets_bed]
    }

    // Run HAPPY for all VCF compared to GIAB truth
    HAPPY_HAPPY_single(ch_vcf_giab, ch_fasta, ch_fasta_fai, empty, empty, empty)

    // Retrieve true-positives from pairwise comparisons.
    HAPPY_HAPPY_pairwise(ch_vcf_pairwise, ch_fasta, ch_fasta_fai, empty, empty, empty)
    ch_pairwise_vcf_index = HAPPY_HAPPY_pairwise.out.vcf
        .map(addTmpId)
        .join(HAPPY_HAPPY_pairwise.out.tbi.map(addTmpId), by: 0)
        .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, []]}

    // Remove nocall  on VCF + index
    GATK4_SELECTVARIANTS_NOCALL(
        HAPPY_HAPPY_pairwise.out.vcf
            .map(addTmpId)
            .join(HAPPY_HAPPY_pairwise.out.tbi.map(addTmpId), by: 0)
            .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, []]}
    )

    // SelectVariants on VCF + index to select true-positives
    GATK4_SELECTVARIANTS_TP(
        GATK4_SELECTVARIANTS_NOCALL.out.vcf
            .map(addTmpId)
            .join(GATK4_SELECTVARIANTS_NOCALL.out.tbi.map(addTmpId), by: 0)
            .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, []]}
    )

    /*
    BCFTOOLS FILTER to remove filter status from pairwise VCF as overlapping
    variants between two VCFs could be regarded as high confident.
      - removing filter status results in similar results for A-B and B-A comparisons.
    */

    BCFTOOLS_ANNOTATE(
        GATK4_SELECTVARIANTS_TP.out.vcf
            .map(addTmpId)
            .join(GATK4_SELECTVARIANTS_TP.out.tbi.map(addTmpId), by: 0)
            .map{id, meta_vcf, vcf, meta_index, index -> [meta_vcf, vcf, index, [], []]}
        , Channel.empty().toList()
    )

    // Run HAPPY on pairwise true-positives against GIAB truth
    HAPPY_HAPPY_tp_giab(
        BCFTOOLS_ANNOTATE.out.vcf.combine(BCFTOOLS_NORM_GIAB.out.vcf).map(createHappyInput),
        ch_fasta, ch_fasta_fai, empty, empty, empty
    )   

    EditSummaryFileHappy(
        Channel.empty().mix(
            HAPPY_HAPPY_single.out.summary_csv,
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
            EditSummaryFileHappy.out.indel_all_csv,
            EditSummaryFileHappy.out.snp_all_csv,
            EditSummaryFileHappy.out.indel_pass_csv,
            EditSummaryFileHappy.out.snp_pass_csv,
            CheckQC.out.qc_output,
            Workflow_ExportParams.out
        ).collect(), 
        multiqc_yaml, [], []
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "GIABEval Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/multiqc_report.html")
    } else {
        def subject = "GIABEval Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
