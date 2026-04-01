#!/usr/bin/env nextflow
// Include processes, alphabetic order of process alias
include { BCFTOOLS_ANNOTATE } from './modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_INPUT } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_GIAB } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_PAIRWISE_TP } from './modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_PAIRWISE } from './modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SINGLE } from './modules/nf-core/bcftools/reheader/main'
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



// Add banner via log.info
// log.info """\
    // G I A B E V A L   P I P E L I N E
    // ===================================
    // input: ${params.vcf_path}
    // output: ${params.outdir}
    // nist version: ${params.nist_version_to_use}
    // assembly: ${params.assembly[params[params.nist_version_to_use].assembly]}
    // GIAB settings: ${params[params.nist_version_to_use]}
    // ===================================

    // """
// .stripIndent(true)

workflow {
    def createMetaWithIdName = {file -> [[id: file.getSimpleName()], file]}
    def addTmpId = {meta, file -> [meta.id, meta, file]}
    def createHappyInput = {meta_query, query, meta_truth, truth ->
        def meta = [
            id: meta_query.id + "_" + meta_truth.id,
            query: meta_query.id,
            truth: meta_truth.id
        ]
        return [meta, query, truth, regions_bed, []]
    }

    def analysis_id = params.outdir.split('/')[-1]

    genome_config = params.assembly[params.genome_build]
    truthset_config = params.truthsets[params.genome_build][params.nist_version]

    // Reference file channels

    ch_fasta = channel.fromPath(genome_config.ref_fasta)
        .map(createMetaWithIdName)
        .first()
    ch_fai = channel.fromPath(genome_config.ref_fai)
        .map(createMetaWithIdName)
        .first()

    ch_false_positives_bed = channel.fromPath(truthset_config.false_positives_bed)
        .map(createMetaWithIdName)
        .first()

    // Reference bed files
    regions_bed = genome_config.exome_target_bed

    // GIAB reference file channels
    ch_giab_truth = Channel.fromPath(truthset_config.truth_vcf)
    .map{file ->
        def samplename = file.name.tokenize("_")
            [[id: "${samplename[0]}_truth"], file]
    }
    .first()

    // Input vcf file channel
    ch_vcf_files = Channel.fromPath(["${params.vcf_path}/*.vcf.gz", "${params.vcf_path}/*.vcf"])
    .map { vcf ->
        // Split filename using params.delim and select indices to create unique identifier
        def tokens = vcf.name.tokenize(params.delim)
        def id_items = params.id_index.collect{idx -> tokens[idx]}
        def identifier = (id_items.join("_") ?: id_items)
        def meta = [
            id: identifier,
            vcf: vcf.simpleName,
            single_end:false
        ]
        [meta, vcf]
    }

    //Compress all VCF and index
    BCFTOOLS_VIEW_INPUT(
        ch_vcf_files
            .map { meta, vcf -> [ meta, vcf, [] ] },
        [],
        [],
        []
    )
    BCFTOOLS_VIEW_GIAB(
        ch_giab_truth
            .map { meta, vcf -> [ meta, vcf, [] ] },
        [],
        [],
        []
    )

    // Slice Input VCFs for primary contigs
    BCFTOOLS_VIEW_PRIMARY(
        BCFTOOLS_VIEW_INPUT.out.vcf
            .join(BCFTOOLS_VIEW_INPUT.out.tbi)
            .map { meta, vcf, tbi -> [ meta, vcf, tbi ]},
        [],
        [],
        []
    )

    /*
    BCFTOOLS_NORM (normalisation) is required to
        - place an indel at the left-most position (left-align)
        - normalizes split multiallelic sites into biallelics
    */
    BCFTOOLS_NORM_INPUT(
        BCFTOOLS_VIEW_PRIMARY.out.vcf
            .join(BCFTOOLS_VIEW_PRIMARY.out.tbi),
        ch_fasta
    )
    BCFTOOLS_NORM_GIAB(
        ch_giab_truth
            .join(BCFTOOLS_VIEW_GIAB.out.tbi),
        ch_fasta
    )

    // Create a channel with vcfs against giab.
    ch_vcf_giab = BCFTOOLS_NORM_INPUT.out.vcf
        .combine(BCFTOOLS_NORM_GIAB.out.vcf)
        .map(createHappyInput)

    // Run HAPPY for all VCF compared to GIAB truth
    HAPPY_HAPPY_single(
        ch_vcf_giab,
        ch_fasta,
        ch_fai,
        ch_false_positives_bed,
        [[:],[]],
        [[:],[]]
    )

    // Reheader Happy output VCF with reference genome .fai
    BCFTOOLS_REHEADER_SINGLE(
        HAPPY_HAPPY_single.out.vcf.map{meta, vcf -> [meta, vcf, [], []]},
        ch_fasta
            .join(ch_fai)
    )

    if (params.run_pairwise) {

        ch_vcf_pairwise = BCFTOOLS_NORM_INPUT.out.vcf
            .combine(BCFTOOLS_NORM_INPUT.out.vcf)
            .filter { meta1, vcf1, meta2, vcf2 ->
                meta1.id < meta2.id
            }
            .map { meta1, vcf1, meta2, vcf2 ->
                def meta = [
                    id:    "pairwise_${meta1.id}_${meta2.id}",
                    query: meta1.id,
                    truth: meta2.id
                ]
                [meta, vcf1, vcf2, regions_bed, []]
            }


        // Retrieve true-positives from pairwise comparisons.
        HAPPY_HAPPY_pairwise(
            ch_vcf_pairwise,
            ch_fasta,
            ch_fai,
            ch_false_positives_bed,
            [[:],[]],
            [[:],[]]
        )

        // ch_pairwise_vcf_index = HAPPY_HAPPY_pairwise.out.vcf
            // .map(addTmpId)
            // .join(HAPPY_HAPPY_pairwise.out.tbi.map(addTmpId), by: 0)
            // .map { _id, meta_vcf, vcf, _meta_index, index -> [meta_vcf, vcf, index, []] }

        // Reheader VCF with reference genome .fai
        BCFTOOLS_REHEADER_PAIRWISE_TP(
            HAPPY_HAPPY_pairwise.out.vcf
                .map{ meta, vcf -> [meta, vcf, [], []] },
            ch_fasta
                .join(ch_fai)
        )

        // Remove nocall  on VCF + index
        GATK4_SELECTVARIANTS_NOCALL(
            BCFTOOLS_REHEADER_PAIRWISE_TP.out.vcf
                .map(addTmpId)
                .join(BCFTOOLS_REHEADER_PAIRWISE_TP.out.index
                      .map(addTmpId),
                      by: 0 )
                .map{ _id, meta_vcf, vcf, _meta_index, index -> [meta_vcf, vcf, index, []] }
        )

        // SelectVariants on VCF + index to select true-positives
        GATK4_SELECTVARIANTS_TP(
           GATK4_SELECTVARIANTS_NOCALL.out.vcf
                .map(addTmpId)
                .join(GATK4_SELECTVARIANTS_NOCALL.out.tbi
                        .map(addTmpId),
                      by: 0)
                .map{ _id, meta_vcf, vcf, _meta_index, index -> [meta_vcf, vcf, index, []] }
        )

        /*
        TODO: BCFTOOLS FILTER to remove filter status from pairwise VCF
        overlapping variants between two VCFs could be regarded as high confident
        thus might not need a filter status (/PASS for all)
          - removing filter status results in similar results for A-B and B-A comparisons
            which is currently not the case.
        */

        BCFTOOLS_ANNOTATE(
            GATK4_SELECTVARIANTS_TP.out.vcf
                .map(addTmpId)
                .join(GATK4_SELECTVARIANTS_TP.out.tbi
                      .map(addTmpId),
                      by: 0)
                .map{ _id, meta_vcf, vcf, _meta_index, index -> [meta_vcf, vcf, index, [], []] },
            []
        )

        // Run HAPPY on pairwise true-positives against GIAB truth
        HAPPY_HAPPY_tp_giab(
            BCFTOOLS_ANNOTATE.out.vcf
                .combine(BCFTOOLS_NORM_GIAB.out.vcf)
                .map(createHappyInput),
            ch_fasta,
            ch_fai,
            ch_false_positives_bed,
            [[:],[]],
            [[:],[]]
        )

        // Reheader Happy pairwise VCF with reference genome .fai
        BCFTOOLS_REHEADER_PAIRWISE(
            HAPPY_HAPPY_tp_giab.out.vcf.map{meta, vcf -> [meta, vcf, [], []]},
            ch_fasta
                .join(ch_fai)
        )

    }

    // Make empty channel if run_pairwise is false
    pairwise_summary = params.run_pairwise ? HAPPY_HAPPY_tp_giab.out.summary_csv : []

    EditSummaryFileHappy(
        Channel.empty().mix(
            HAPPY_HAPPY_single.out.summary_csv,
            pairwise_summary
        ).ifEmpty([])

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
