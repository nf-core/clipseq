process GET_CONSENSUS_COUNTS {
    tag "${meta.id}"
    label "process_high"

    conda "bioconda::bioconductor-rtracklayer=1.62.0"
    container 'quay.io/biocontainers/bioconductor-rtracklayer:1.62.0--r43ha9d7317_0'

    input:
    tuple val(meta), path(bedtools_map_outputs)
    tuple val(meta2), path(gtf)
    val(output_name)

    output:
    tuple val(meta), path ("*.tsv")         , emit: tsv
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    template 'prepare_consensus_count_table.R'
}
