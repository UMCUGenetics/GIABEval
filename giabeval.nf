#!/usr/bin/env nextflow
include { RTGTOOLS_VCFEVAL } from './modules/modules/nf-core/rtgtools/vcfeval/main.nf'
include { HAPPY_HAPPY } from "./modules/modules/nf-core/happy/happy/main.nf"

workflow {
    rtg_input = Channel.of([[id: "test"],
            "/hpc/diaggen/projects/GIAB_eval/Test/SNP_nistregions_230623_A00295_0741_AH2LJWDSX7_Val_Crev4_GiabU175754_JAN8898.vcf.gz",
            "/hpc/diaggen/projects/GIAB_eval/Test/SNP_nistregions_230623_A00295_0741_AH2LJWDSX7_Val_Crev4_GiabU175754_JAN8898.vcf.gz.tbi",
            "/hpc/diaggen/projects/GIAB_eval/GIAB_database/hg19/HG001/v2.19/SNV_NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz",
            "/hpc/diaggen/projects/GIAB_eval/GIAB_database/hg19/HG001/v2.19/SNV_NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz.tbi",
            "/hpc/diaggen/projects/CREv4_analysis/Dx_tracks/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank.bed",
            []
    ])
    rtg_sdf = Channel.of([[id: "test_sdf"], "/hpc/diaggen/projects/GIAB_eval/ref_genome_index/Homo_sapiens.GRCh37.GATK.illumina_SDF"])
    RTGTOOLS_VCFEVAL(rtg_input, rtg_sdf)

}

workflow.onComplete {}



