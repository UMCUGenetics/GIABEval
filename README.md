

[![GitHub Actions CI Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20CI/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20linting/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# GIABEval
Genome-in-a-bottle evaluation software to determine statistics such as precision, recall and f-measure.

Note that statistics calculations will be on the primary contigs as defined in nextflow.config.
This prevents issues due to diffent reference genomes such als ALT, Decoy, or tertairy software reference genomes.

## Installation

An install script (`install.sh`) is available, tailored to the UMCU infrastructure. For other users, we recommend following these steps:

1. **Install Nextflow** by following the official guide: [https://docs.seqera.io/nextflow/install](https://docs.seqera.io/nextflow/install)

2. **Clone the GIABEval repository** and initialise the submodules:
```bash
   git clone git@github.com:UMCUGenetics/GIABEval.git
   cd GIABEval
   git submodule update --init --recursive
```


## Data files

GIABEval requires a few specific reference input files.

- Genome reference files (GRCh38): https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
- The `rtg_index` may be downloaded from [realtimegenomics](https://www.realtimegenomics.com/news/pre-formatted-reference-datasets) or can alternatively produced using RTG Tools with command: `rtg format -o reference.sdf reference.fasta`, using the genome fasta as input
- The `exome_target_bed` is available for [GRCh37](https://github.com/UMCUGenetics/Dx_tracks/blob/master/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank_collapsed.bed) and [GRCh38](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Functional/GRCh38_refseq_cds.bed.gz).

The truthsets can be obtained from the ncbi, for example HG002: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/


## Configuration

GIABEval requires reference files for the genome assembly and GIAB truthsets. By default, these point to locations on the UMCU cluster, as defined in `conf/genomes_truthsets.config`. To run GIABEval elsewhere, you'll need to override these paths with your own custom config file.

### Creating a custom config

Create a file (e.g. `my_resources.config`) following the structure below. You only need to include the assemblies and truthsets you want to use — anything you omit will fall back to the defaults.
```groovy
params {
  assembly {
    "GRCh38" {                              // assembly name — referenced by truthsets below
      ref_fasta        = "/path/to/genome.fna"
      ref_fai          = "/path/to/genome.fna.fai"
      rtg_index        = "/path/to/genome.SDF/"
      exome_target_bed = "/path/to/targets.bed"
      primary_contigs  = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    }
  }

  truthsets {
    "GRCh38" {                              // must match an assembly name above
      "hg002_nist_v4_2_1" {                 // truthset ID — pick any name you like
        input_base          = "/path/to/GIAB/HG002/NISTv4.2.1/GRCh38"
        truth_vcf           = "${input_base}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        truth_vcf_index     = "${input_base}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
        false_positives_bed = "${input_base}/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
      }
    }
  }
}
```
Then pass it to Nextflow with `-c my_resources.config` alongside with the required parameters `--genome_build GRCh38` `--nist_version hg002_nist_v4_2_1` when running the pipeline.

### Nextflow profiles

By default, GIABEval submits jobs via **SLURM**. When running outside the UMCU cluster, you'll need to adjust the institute-specific SLURM settings — most notably:

- `cluster_account`
- `singularity_cachedir`
- `singularity_runoptions`

To run GIABEval locally or on a non-SLURM scheduler, you can define your own Nextflow profile. See the Nextflow documentation for details:

- [Config profiles](https://docs.seqera.io/nextflow/config#config-profiles)
- [Executors](https://docs.seqera.io/nextflow/executor)


## Usage
### UMCU

```bash
nextflow run GIABeval/main.nf \
  --vcf_path <input_vcf_dir_path/> \
  --outdir <output_dir_path> \
  --email <email>
```

### Other institute SLURM

``` bash
nextflow run GIABeval/main.nf \
  -c my_resources.config \
  --genome_build GRCh38> \ # modify accordingly
  --nist_version hg002_nist_v4_2_1 \ # modify accordingly
  --cluster_account <account_name> \
  --singularity_cachedir </path/to/singularity/cachedir> \
  --singularity_runoptions "" \
  --vcf_path <input_vcf_dir_path/> \
  --outdir <output_dir_path> \
  --email <email>
```





## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  GIABEval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
22
