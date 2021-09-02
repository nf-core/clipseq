# nf-core/clipseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Enhancements & fixes

## 1.1.0dev

* Bug fix for RSeQC input channel if skipping deduplication
* Changed parameter `--deduplicate` (default `true`) to `--skip_deduplication` (default `false`)
* Added error strategy as `ignore` for Preseq (need to update to a `--skip_preseq` parameter)
* Fixed `Cutadapt` bug if `sample_id` and `fastq` basename are the same
* Updated Template to nf-core/tools 1.14

## [1.0.0] - 2021-04-27

Initial release of nf-core/clipseq, created with the [nf-core](https://nf-co.re/) template.

### Pipeline summary

1. Adapter and quality trimming (`Cutadapt`)
2. Pre-mapping to e.g. rRNA and tRNA sequences (`Bowtie 2`)
3. Genome mapping (`STAR`)
4. UMI-based deduplication (`UMI-tools`)
5. Crosslink identification (`BEDTools`)
6. Bedgraph coverage track generation (`BEDTools`)
7. Peak calling (multiple options):
    * `iCount`
    * `Paraclu`
    * `PureCLIP`
    * `Piranha`
8. Motif detection (`DREME`)
9. Quality control:
    * Sequencing quality control (`FastQC`)
    * Library complexity (`Preseq`)
    * Regional distribution (`RSeQC`)
10. Overall pipeline run and QC summaries and peak calling comparisons (`MultiQC`)
