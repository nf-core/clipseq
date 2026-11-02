process IDENTIFY_UNBINNED {
    tag "$meta.id"
    label 'process_low'
    label 'process_high_memory'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dde1efad90c922a0198cae64c642be95605c23cf1e53e3c35491817bf6c48b/data':
      'community.wave.seqera.io/library/bedtools_pybedtools_pysam_matplotlib_pruned:79786472f5bd377f' }"

    input:
    tuple val(meta) , path(bam)
    path(discarded)

    output:
    tuple val(meta), path("*discarded_reads_sorted.bam"),  emit: bam
    tuple val(meta), path("*discarded_reads_sorted.bam.bai"),  emit: bai, optional : true
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

     python ${projectDir}/modules/local/identify_unbinned/templates/identify_unassigned.py \
     -b ${bam} \
     -d ${discarded}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

}
