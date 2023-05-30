# nf-core/clipseq: Output

## Introduction

This document describes the output produced by the pipeline. The plots are taken from the MultiQC report, which summarises results at the end of the pipeline and also includes CLIP-specific summary metrics.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the steps described in the main [README.md](https://github.com/nf-core/clipseq/README.md).

* [Preprocessing](#preprocessing)
  * [FastQC](#fastqc) - Raw read QC
  * [UMI-tools](#umi-tools-extract) - Extract UMI from FastQ read and append to read name
  * [Cutadapt](#cutadapt) - Adapter and quality trimming
* [Alignment](#alignment)
  * [Bowtie 2](#bowtie-2) - Pre-alignment to rRNA and tRNA sequences
  * [STAR](#star) - Splice-aware genome alignment
* [Crosslink identification](#crosslink-identification)
  * [UMI-tools](#umi-tools-dedup) - UMI-based PCR deduplication
  * [BEDTools](#bedtools) - Single-nucleotide crosslink position identification
* [Peak calling](#peak-calling)
  * [iCount](#icount)
  * [Paraclu](#paraclu)
  * [PureCLIP](#pureclip)
  * [Piranha](#piranha)
* [Motif finding](#motif-finding)
  * [DREME](#dreme)
* [CLIP summary and quality control](#clip-summary-and-quality-control)
  * [CLIP summary](#clip-summary)
  * [Preseq](#preseq) - Library complexity
  * [RSeQC](#resqc) - Crosslink distribution over genomic features
  * [MultiQC](#multiqc) - CLIP summary metrics and QC for raw reads, alignment, PCR deduplication, library complexity and crosslink distribution over genomic features

## Preprocessing

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory:** `fastqc`

* `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
* `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### UMI-tools extract

Optionally, [UMI-tools](https://umi-tools.readthedocs.io/en/latest/) is used to extract the UMI from the beginning of the read and append it to the read name using `_` as a delimiter. This has often already been done as part of demultiplexing an Illumina sequencing run.

> **NB** As each FASTQ contains only one sample, the whole barcode (experimental and UMI) can be treated as the UMI as described [here](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html#step-3-extract-the-umis).

**Output directory:** `umi`

* `sample.umi.fastq.gz`: FASTQ file of reads with UMIs extracted

### Cutadapt

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) removes adapters and also quality trims the data. By default the pipeline trims the Illumina universal adapter sequences and filters out reads that are shorter that 12 nt after trimming

**Output directory:** `cutadapt`

* `sample.trimmed.fastq.gz`: FASTQ file after trimming
* `sample.cutadapt.log`: Cutadapt log file

## Alignment

### Bowtie 2

For CLIP data analysis it is often important to pre-map to rRNA and tRNA sequences. FASTA files for a number of organisms are provided as part of the pipeline. [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is used to identify these reads.

**Output directory:** `premap`

* `sample.premapped.bam`: BAM file of reads mapped to the premapping index
* `sample.premapped.bam.bai`: BAI file for BAM
* `sample.premap.log`: Premapping (Bowtie 2) log file
* `sample.unmapped.fastq.gz`: FASTQ file of reads that do not map to the premapping index that is passed to the next step of the pipeline.

### STAR

[STAR](https://github.com/alexdobin/STAR) is used to align to the genome. Importantly, soft-clipping of the 5' end of the read is prevented, ensuring the crosslink position can be correctly identified.

**Output directory:** `mapped`

* `sample.Aligned.sortedByCoord.bam`: BAM file of reads mapped to the genome
* `sample.Aligned.sortedByCoord.bam.bai`: BAI file for BAM
* `sample.Log.final.out`: Alignment (STAR) log file

## Crosslink identification

### UMI-tools dedup

[UMI-tools](https://umi-tools.readthedocs.io/en/latest/) is used for UMI aware PCR deduplication. The directional method is used.

**Output directory:** `dedup`

* `sample.dedup.bam`: BAM file of deduplicated reads
* `sample.dedup.bam.bai`: BAI file for BAM
* `sample.log`: Deduplication (UMI-tools) log file

### BEDTools

[BEDTools](https://bedtools.readthedocs.io/en/latest/) is used to identify the crosslinks from the BAM files. The crosslink BED files are single-nucleotide resolution (i.e. each entry is 1 nt wide) and the score is the number of crosslinks at that position. In the crosslink BEDGRAPH files, a positive score indicates a crosslink on the positive strand and a negative score one on the negative strand.

**Output directory:** `xlinks`

* `sample.xl.bed.gz`: BED file of crosslinks
* `sample.xl.bedgraph.gz`: BEDGRAPH file of crosslinks

## Peak calling

The following peak callers are currently provided in the pipeline:

* [iCount](https://icount.readthedocs.io/en/latest/)
* [Paraclu](http://cbrc3.cbrc.jp/~martin/paraclu/)
* [PureCLIP](https://pureclip.readthedocs.io/en/latest/)
* [Piranha](https://github.com/smithlabcode/piranha)

The user can specify which one(s) are run. Filenames with the default run parameters are shown below, but are adjusted by the pipeline according to the parameters specified.

### iCount

**Output directory** `icount`

* `sample.3nt.sigxl.bed.gz`: BED file of significant crosslink positions using a 3 nt half-window setting
* `sample.3nt_3nt.peaks.bed.gz` BED file of peaks using a 3 nt half window and a 3 nt merge window

### Paraclu

**Output directory** `paraclu`

* `sample.10_200nt_2.peaks.bed.gz`: BED file of peaks using a minimum value/score of 10, a maximum cluster length of 200 and a minimum density increase of 2.

### PureCLIP

**Output directory** `pureclip`

* `sample.sigxl.bed.gz`: BED file of significant crosslink sites
* `sample.8nt.peaks.bed.gz`: BED file of peaks using a merge distance of 8 nt

### Piranha

**Output directory** `piranha`

* `sample.3nt_3nt.peaks.bed.gz`: BED file of peaks using a bin size of 3 and a cluster distance of 3

## Motif finding

### DREME

[DREME](http://meme-suite.org/doc/dreme.html) is used for basic motif calling used peaks. By default the sequence from the region +/- 20 nt around the crosslink site is provided as input for DREME.

**Output directories** `icount_motif`, `paraclu_motif`, `pureclip_motif`, `piranha_motif`

* `sample_dreme/`: Directory containing DREME output files: `dreme.html`, `dreme.txt`, `dreme.xml`

## CLIP summary and quality control

### CLIP summary

The pipeline uses a custom script to provide CLIP-specific summary metrics that are plotted in the MultiQC report.

#### Mapping

This section plots the counts/percentages of reads mapped to the premapping index, mapped to the genome, and that remain unmapped.

#### Deduplication

This section plots three measures from the UMI-based PCR deduplication.

  1. *Reads* shows the number of reads before and after deduplication.
  2. *Ratios* shows the PCR deduplication ratio.
  3. *Mean UMIs* shows the mean number of unique UMIs per position.

#### Crosslinks

This section plots two measures from crosslink identification.

  1. *Counts* shows the number of crosslinks and crosslink sites.
  2. *Ratios* shows the ratio of crosslinks to crosslink sites.

#### Peaks

This sections plots three peak-calling metrics (if peak calling has been performed) to enable comparison of different tools and optimisation of specific peak-caller parameters.

  1. *Crosslinks in peaks* shows the total percentage of crosslinks within peaks.
  2. *Crosslink sites in peaks* shows the total percentage of crosslink sites within peaks.
  3. *Peak-crosslink coverage* shows the total percentage of nucleotides within peaks that are covered by a crosslink site.

**Output directory** `clipqc`

* `*.tsv`: TSV files containing derived metrics from the pipeline outputs used to produce the MultiQC CLIP summary metric plots.

### Preseq

[Preseq](http://smithlabresearch.org/software/preseq/) is used to estimate the complexity of the sequenced library.

**Output directory** `preseq`

* `sample.ccurve.txt`: TXT file of complexity curve Preseq output.
* `sample.command.log`: Preseq log file

### RSeQC

[RSeQC](http://rseqc.sourceforge.net/) is used to calculate how deduplicated mapped reads are distributed over genomic features.

**Output directory** `rseqc`

* `sample.read_distribution.txt`: TXT file of `read_distribution.py` output from RSeQC.

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
