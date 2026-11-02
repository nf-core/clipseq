process ENCODE_MOVEUMI {
    label "process_single"

    conda "bioconda::biopython=1.78 pigz=2.6"
    container "quay.io/biocontainers/mulled-v2-877c4e5a8fad685ea5bde487e04924ac447923b9:b7daa641364165419b9a87d9988bc803f913c5b6-0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    process_name = task.process
    template 'encode_moveumi.py'
}
