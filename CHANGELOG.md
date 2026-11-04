# nf-core/clipseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/clipseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

- Replace legacy Groovy template helpers with shared nf-core utilities and pipeline-specific initialization/reporting subworkflows.
- Migrate parameter validation, help, and summaries to JSON Schema draft 2020-12 through nf-schema.
- Configure resource caps with `process.resourceLimits` instead of `max_*` parameters.

### `Dependencies`

- Require Nextflow >=25.10.4 and pin nf-schema 2.7.3.
- Update nf-core tools to 4.2.1 and CI nf-test to 0.9.5.

### `Deprecated`
