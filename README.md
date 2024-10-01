

[![GitHub Actions CI Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20CI/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20linting/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# GIABEval
Genome-in-a-bottle evaluation software to determine statistics such as precision, recall and f-measure.

Note that statistics calculations will be on the primary contigs as defined in nextflow.config.
This prevents issues due to diffent reference genomes such als ALT, Decoy, or tertairy software reference genomes.

## Get submodules, such as NextflowModules, CustomModules and install OpenJDK and Nextflow
```bash
sh install.sh
```

## Usage

```bash
nextflow run main.nf -c nextflow.config --vcf_path [input_vcf_dir_path] --outdir [output_dir_path] --email [email]
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
