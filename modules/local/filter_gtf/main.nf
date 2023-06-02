process CLIPSEQ_FILTER_GTF {
    tag "$gtf"
    label "process_single"

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3':
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: "${gtf.simpleName}_filtered.gtf"
    template 'filter_gtf.py'
}
