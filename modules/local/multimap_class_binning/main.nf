process MULTIMAP_CLASS_BINNING {
    label "process_medium"
    label "process_high_memory"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dde1efad90c922a0198cae64c642be95605c23cf1e53e3c35491817bf6c48b/data':
    'community.wave.seqera.io/library/bedtools_pybedtools_pysam_matplotlib_pruned:79786472f5bd377f' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta_gtf), path(gtf)

    output:
    path "discarded_ambiguous_reads*.tsv", emit: discarded_reads
    path("*_intersected.bed"), emit: bed
    tuple val (meta), path("repeat_ncRNA_read_counts_*.csv"), emit: all_reads
    tuple val (meta), path("repeat_ncRNA_read_counts5*.csv"), emit: counted_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    python ${projectDir}/modules/local/multimap_class_binning/templates/multimap_class_binning.py \
    --bamdir . \
    --gtf ${gtf} \
    -i ${bam} \
    -o . \
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pandas: \$(python -c "import pybedtools; print(pybedtools.__version__)")
        pandas: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pandas: \$(python -c "import pybedtools; print(pybedtools.__version__)")
        pandas: \$(python -c "import pysam; print(pysam.__version__)")
        pandas: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

}
