// modules/local/bedgraph_strand_split/main.nf
process BEDGRAPH_STRAND_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bedgraph)

    output:
    tuple val(meta), path("*.pos.bedgraph"), emit: pos_bedgraph
    tuple val(meta), path("*.neg.bedgraph"), emit: neg_bedgraph
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract positive and negative strands (4th column > 0 and < 0)
    awk -v OFS='\t' '\$4 > 0 {print \$1, \$2, \$3, \$4}' ${bedgraph} > ${prefix}.pos.bedgraph
    awk -v OFS='\t' '\$4 < 0 {print \$1, \$2, \$3, -\$4}' ${bedgraph} > ${prefix}.neg.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pos.bedgraph
    touch ${prefix}.neg.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}