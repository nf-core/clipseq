# ![nf-core/clipseq](docs/images/nf-core-clipseq_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/clipseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/clipseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/clipseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/clipseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/clipseq.svg)](https://hub.docker.com/r/nfcore/clipseq)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23clipseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/clipseq)

## Introduction

**nf-core/clipseq** is a bioinformatics best-practice analysis pipeline for CLIP (cross-linking and immunoprecipitation) sequencing data analysis to study RNA-protein interactions.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary

By default, the pipeline currently performs the following:

1. Adapter and quality trimming (`Cutadapt`)
2. Pre-mapping to e.g. rRNA and tRNA sequences (`Bowtie 2`)
3. Genome mapping (`STAR`)
4. UMI-based deduplication (`UMI-tools`)
5. Crosslink identification (`BEDTools`)
6. Bedgraph coverage track generation (`BEDTools`)
7. Peak calling (multiple options):
    - `iCount`
    - `Paraclu`
    - `PureCLIP`
    - `Piranha`
8. Motif detection (`DREME`)
9. Quality control:
    - Sequencing quality control (`FastQC`)
    - Library complexity (`Preseq`)
    - Regional distribution (`RSeQC`)
10. Overall pipeline run and QC summaries and peak calling comparisons (`MultiQC`)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/clipseq -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/clipseq -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '[path to design file]' --fasta '[path to genome FASTA]'
    ```

See [usage docs](https://nf-co.re/clipseq/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/clipseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/clipseq/usage) and [output](https://nf-co.re/clipseq/output).

## Credits

nf-core/clipseq was originally written by Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) and Anob Chakrabarti ([@amchakra](https://github.com/amchakra)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK.

It started life in April 2020 as a Nextflow DSL2 Luscombe Lab ([@luslab](https://github.com/luslab)) lockdown hackathon day and we thank all the lab members for their early contributions.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#clipseq` channel](https://nfcore.slack.com/channels/clipseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/clipseq for your analysis, please cite it using the following doi: [10.5281/zenodo.4723016](https://doi.org/10.5281/zenodo.4723016)

References of tools and data used in this pipeline can be found in [CITATIONS.md](https://github.com/nf-core/clipseq/CITATIONS.md)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
