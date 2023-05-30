# nf-core/clipseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/clipseq/usage](https://nf-co.re/clipseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Main arguments

### `--input`

You will need to create a design file with information about the samples in your experiment before running the pipeline. Only single end reads are currently supported. Use this parameter to specify its location.

```bash
--input '[path to design file]'
```

It has to be a comma-separated file with 2 columns, and a header row as shown in the examples below. The column headers must be `sample` and `fastq`. By naming the `sample` rows uniquely, one can identify and simultaneously run multiple replicates and samples:

```bash
sample,fastq
exp1_rep1,clip0001_01.fastq.gz
exp1_rep2,clip0001_02.fastq.gz
exp2_rep1,clip0002_01.fastq.gz
exp2_rep2,clip0002_02.fastq.gz
```

| Column         | Description                                                                                                 |
|----------------|-------------------------------------------------------------------------------------------------------------|
| `sample`        | Unique identifier for read, which may include information about sample and replicate. |
| `fastq`    | Full path to FastQ file for read. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz". |

### `--fasta`

Full path to fasta file containing reference genome (mandatory if --genome is not specified). If you don't have a STAR index available this will be generated for you automatically. Alternatively, it can be set using `--star_index`.

```bash
--fasta '[path to FASTA reference]'
```

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag. If you have the iGenomes references locally available you can set `--igenome_base`, otherwise they will be automatically obtained from AWS-iGenomes. You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* _H. sapiens_
  * `--genome GRCh37`
* _M. musculus_
  * `--genome GRCm38`
* _D. melanogaster_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

Premapping to rRNA and tRNA will be automatically triggered if there is a reference available for the iGenomes reference chosen. See [smallRNA config file](../conf/smrna.config) for availability. Alternatively, this can be set by `--smrna_org` or `--smrna_fasta` as shown below.

### `--smrna_org`

The pipeline comes equipped with some 'smallRNA' FASTA references for premapping. This includes rRNA and tRNA sequences, the sources of which can be viewed [here](https://github.com/ulelab/smallRNAs/blob/main/Sources.md). The purpose of this premapping is to capture abundant ncRNA that are present in multiple similar copies in the genome, making them hard to assign reads to. tRNA can occur within genes and without proper handling can result in misassignment of reads to mRNA in [certain situations](https://rnajournal.cshlp.org/content/early/2018/08/21/rna.067348.118). The user may also be interested in the tRNA and rRNA binding of their protein, and this premapping enables simple assessment of this binding. These are available for the following organisms:

* Human
  * `--smrna_org human`
* Mouse
  * `--smrna_org mouse`
* Rat
  * `--smrna_org rat`
* Zebrafish
  * `--smrna_org zebrafish`
* Fruitfly
  * `--smrna_org fruitfly`
* Yeast
  * `--smrna_org yest`

### `--smrna_fasta`

Alternatively, the RNA premapping reference can be supplied by the user by giving the path to the reference FASTA:

```bash
--smrna_fasta '[path to smallRNA FASTA reference]'
```

### `--peakcaller`

By default, peak calling on identified crosslinks is not performed unless specified using the `--peakcaller` argument. Currently the following peak callers are implemented:

* [iCount](https://icount.readthedocs.io/en/latest/)
  * `--peakcaller icount`
* [Paraclu](http://cbrc3.cbrc.jp/~martin/paraclu/)
  * `--peakcaller paraclu`
* [PureCLIP](https://pureclip.readthedocs.io/en/latest/)
  * `--peakcaller pureclip`
* [Piranha](https://github.com/smithlabcode/piranha)
  * `--peakcaller piranha`

Multiple peak callers can specified separated using a comma, e.g. `--peakcaller icount,paraclu`. As a short-hand, all peak callers can be specified using `--peakcaller all`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/clipseq --input '[path to design file]' --fasta '[path to genome FASTA]' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/clipseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/clipseq releases page](https://github.com/nf-core/clipseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/clipseq`](https://hub.docker.com/r/nfcore/clipseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/clipseq`](https://hub.docker.com/r/nfcore/clipseq/)
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/clipseq`](https://hub.docker.com/r/nfcore/clipseq/)
* `shifter`
  * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
  * Pulls software from Docker Hub: [`nfcore/clipseq`](https://hub.docker.com/r/nfcore/clipseq/)
* `charliecloud`
  * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
  * Pulls software from Docker Hub: [`nfcore/clipseq`](https://hub.docker.com/r/nfcore/clipseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
