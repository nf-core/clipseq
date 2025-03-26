<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-clipseq_logo_dark.png">
    <img alt="nf-core/clipseq" src="docs/images/nf-core-clipseq_logo_light.png">
  </picture>
</h1>

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/clipseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![GitHub Actions CI Status](https://github.com/nf-core/clipseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/clipseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/clipseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/clipseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/clipseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/clipseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23clipseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/clipseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/clipseq** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

By default, the pipeline currently performs the following:

1. Adapter and quality trimming (`TrimGalore!`)
2. Pre-mapping to e.g. rRNA and tRNA sequences (`Bowtie`)
3. Genome mapping (`STAR`)
4. UMI-based deduplication (`UMI-tools`)
5. Crosslink identification (`BEDTools`)
6. Bedgraph coverage track generation (`BEDTools`)
7. Crosslink summaries over RNA types (`iCount Summary`)
8. Crosslink metagene profiles (`iCount Metagene`)
9. Peak calling (multiple options):
   - `iCount`
   - `Paraclu`
   - `PureCLIP`
   - `Clippy`
10. Consensus peak output, ie. grouping all crosslinks, calling peaks and reporting a table with individual sample counts over consensus peaks. Useful as an input for differential binding analyses.
11. Motif detection (`PEKA`)
12. Overall pipeline run and QC summaries and peak calling comparisons (`MultiQC`)

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample_name,group_name,input_name,fastq
PHO92_A,PHO92,HNRNPC,https://github.com/luslab/test-datasets/raw/clipseq/v_2_0/fastq/ERR3988069-yeast-quarter.fastq.gz
PHO92_B,PHO92,HNRNPC,https://github.com/luslab/test-datasets/raw/clipseq/v_2_0/fastq/ERR3988069-yeast-quarter.fastq.gz
PHO92_C,PHO92,HNRNPC,https://github.com/luslab/test-datasets/raw/clipseq/v_2_0/fastq/ERR3988069-yeast-quarter.fastq.gz
HNRNPC,,,https://github.com/luslab/test-datasets/raw/clipseq/v_2_0/fastq/ERR3988069-yeast-quarter.fastq.gz
```

Each row represents a fastq file (single-end). If multiple rows have the same `sample_name`, `group_name` and `input_name` then their fastqs will be concatenated at the beginning of the pipeline (eg. in the case of a sample that was resequenced and so has multiple fastq files). If multiple rows have different `sample_name`s but the same `group_name` and `input_name`, then the crosslinks for these samples will be merged and peaks and downstream analyses will be run on the grouped crosslinks.

The `input_name` is used for providing input data, currently the only peak caller we support that can use input data is PureCLIP. The `input_name` must match either a `sample_name` or a `group_name`.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/clipseq \
   -profile test, <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://nf-co.re/clipseq/usage) and the [parameter documentation](https://nf-co.re/clipseq/parameters).

## Pipeline output

To see the the results of a test run with a full size dataset refer to the [results](https://nf-co.re/clipseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/clipseq/output).

## A note on paired-end reads

The pipeline currently does not support paired-end reads, as in our experience alignment using both reads when available doesn't improve analysis of CLIP data. When recieving CLIP data sequenced paired-end, we recommend running the pipeline with the read containing the crosslink and ensuring the crosslink_position parameter is set appropriately. If you have evidence to the contrary please do get in touch and let us know, or if you are working on a new variant protocol where paired-end alignment is important please do reach out.

## Credits

nf-core/clipseq was originally written by Charlotte West ([@charlotte-west](https://github.com/charlotte-west)) and Anob Chakrabarti ([@amchakra](https://github.com/amchakra)) from [Luscombe Lab](https://www.crick.ac.uk/research/labs/nicholas-luscombe) at [The Francis Crick Institute](https://www.crick.ac.uk/), London, UK. It started life in April 2020 as a Nextflow DSL2 Luscombe Lab ([@luslab](https://github.com/luslab)) lockdown hackathon day and we thank all the lab members for their early contributions.

v2.0 was spearheaded by [Ule lab](https://github.com/ulelab) in collaboration with [Goodwright Ltd.](https://goodwright.com/); namely Charlotte Capitanchik ([@CharlotteAnne](https://github.com/charlotteanne)), Chris Cheshire and Sam Ireland.

We thank the following people for their extensive assistance in the development of this pipeline:
Ira Iosub, Marc Jones, Rupert Faraway, Oscar G. Wilkins, Klara Kuret, Simon Murray

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#clipseq` channel](https://nfcore.slack.com/channels/clipseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/clipseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
