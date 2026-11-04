---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/clipseq`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/clipseq.

This page describes the recommended nf-core way to contribute to both nf-core/clipseq and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/clipseq.

> [!NOTE]
> If you need help using or modifying nf-core/clipseq, ask on the nf-core Slack [#clipseq](https://nfcore.slack.com/channels/clipseq) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/clipseq/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/clipseq repository](https://github.com/nf-core/clipseq) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focused as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/clipseq documentation](https://nf-co.re/clipseq/usage) or ask on the nf-core [#clipseq](https://nfcore.slack.com/channels/clipseq) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/clipseq without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/clipseq/codespaces).

### Development requirements

Use the following versions to match the repository's development and test configuration:

- **Nextflow 25.10.4**: the minimum version in `nextflow.config` and the minimum tested by CI.
- **Java 17 through 25**: the supported Java range for [Nextflow 25.10.4](https://github.com/nextflow-io/nextflow/blob/v25.10.4/docs/install.md).
- **nf-core/tools 4.2.1**: the version recorded in `.nf-core.yml` and installed by linting CI.
- **nf-test 0.9.5**: the version used by the nf-test CI workflow.
- **Python 3.12**: the version used by linting CI. nf-core/tools 4.2.1 supports Python >=3.10 and <4.

The manifest permits newer Nextflow releases, but check parser compatibility when upgrading.
Nextflow 26.04 and newer enable the stricter v2 parser by default; code using legacy syntax may need `NXF_SYNTAX_PARSER=v1` until migrated.
See the [Nextflow parser-version guidance](https://training.nextflow.io/latest/info/nxf_versions/).

The pipeline pins `nf-schema@2.7.3` in `nextflow.config` for parameter validation.
That plugin requires Nextflow >=25.10.0 and Java >=17, so the versions above meet its requirements.
Install Docker or Apptainer for container-based tests.

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test tests/nf_test/test --profile debug,test,docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Inspect any snapshot differences first. If they are expected consequences of your changes, update the affected snapshots with the following command:

```bash
nf-test test tests/nf_test/test --profile debug,test,docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint --dir .
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/clipseq code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channels if needed. Follow the pipeline-specific version-reporting conventions below and add relevant reports to the MultiQC input channel.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline initialization subworkflow, `subworkflows/local/utils_nfcore_clipseq_pipeline/main.nf`.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `docs/usage.md`, `docs/output.md`, and `CITATIONS.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline-specific contribution guidelines

### CLIPseq development branches

The general nf-core contribution workflow targets `dev`.
For work on the ongoing CLIPseq 2.0 development, coordinate the pull request target with the maintainers; the current integration branch is `feat-2-0`.
Create a focused feature branch from the intended target and keep unrelated changes out of the pull request.

### Workflow and module structure

`main.nf` initializes the pipeline and calls `workflows/clipseq.nf`.
Reusable CLIPseq processing logic belongs in `subworkflows/local/`, and pipeline-specific processes belong in `modules/local/`.
Follow the structure, process aliases, and channel conventions of nearby code.

Use nf-core/tools to install or update shared components under `modules/nf-core/` and `subworkflows/nf-core/`, and include the corresponding `modules.json` changes.
Keep local customization separate from upstream components where possible; record any necessary upstream module patches through nf-core/tools.

### Parameter and samplesheet validation

Define parameter defaults in `nextflow.config` and document their types, defaults, descriptions, and help text in `nextflow_schema.json`.
The parameter schema uses JSON Schema draft 2020-12, with `$defs` and matching references.
Use the pinned `nf-schema` plugin for parameter validation, help, and parameter summaries.
Pipeline-specific initialization checks belong in `subworkflows/local/utils_nfcore_clipseq_pipeline/main.nf`.

Preserve the source-specific samplesheet validation performed by `INPUT_CHECK` and its local samplesheet-checking module.
FASTQ input uses `sample_name,group_name,input_name,fastq` columns; deduplicated BAM input uses `sample_name,group_name,input_name,bam`.
When changing samplesheet handling, test sample merging, replicate groups, control assignments, and validation failures as appropriate.

### Software versions and MultiQC

Updated nf-core modules publish software versions to the `versions` topic.
Local modules that still emit `versions.yml` must have those outputs mixed into `ch_versions`.
The reporting section of `workflows/clipseq.nf` combines both sources for MultiQC.
Do not reintroduce the removed Groovy schema helpers or a separate software-version scraping process.

When adding reports, connect them to the appropriate MultiQC input channels and update `assets/multiqc_config.yml` and `docs/output.md` as needed.
Preserve the workflow summary, methods description, and optional completion reporting.

### Resources and containers

Use the shared resource labels in `conf/base.config` and pass allocated resources to tools through `task.cpus` and `task.memory` where supported.
Resource caps use `process.resourceLimits`; do not add the removed `max_cpus`, `max_memory`, or `max_time` parameters.

Use pinned, fully qualified container references.
The local BioContainers modules use explicit Quay addresses for Docker, Singularity, and Apptainer.
Shared nf-core modules retain their upstream container definitions, including Seqera Containers where supplied.
Verify that the referenced tag exists before changing a container address or version.

### Targeted tests and cluster execution

Existing pipeline tests live in `tests/nf_test/test/`.
Start with tests covering the changed stage, including relevant skip flags and optional workflow branches.
For example, run the preprocessing tests with Docker:

```bash
nf-test test tests/nf_test/test/test_only_preprocessing.nf.test --profile debug,test,docker --verbose
```

Replace `docker` with `apptainer` when testing on a system configured for Apptainer.
For lightweight startup and parameter-validation checks, use:

```bash
nf-test test tests/nf_test/test/schema_validation.nf.test --verbose
```

Run data-processing tests on allocated cluster resources or through CI, rather than on a login node.
When each task should be a separate Slurm job, load a cluster configuration that sets `process.executor = 'slurm'`.
The Nextflow controller and its tasks must be able to access the repository, work directory, inputs, outputs, and container cache.

The `test_full` profile exercises larger scenarios; use it after targeted tests pass and request resources appropriate to the selected data.
Document which tests you ran and any validation that remains incomplete in the pull request.

### Documentation and acknowledgements

Update `CHANGELOG.md` for user-visible changes and dependency updates.
Update `docs/usage.md`, `docs/output.md`, `README.md`, and `CITATIONS.md` where affected.
For contributor acknowledgements, provide the person's preferred name and GitHub handle, and describe their contribution using the existing README conventions.
Record people in `manifest.contributors` in `nextflow.config`, using `name`, `github` where known, and a `contribution` list containing `author`, `maintainer`, or `contributor`.
A person can have multiple roles, for example `contribution: ['author', 'maintainer']`.
Author roles are recorded in the structured contributor entries; keep the `.nf-core.yml` template author list consistent with those entries.
