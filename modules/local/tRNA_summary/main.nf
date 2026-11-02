process tRNA_SUMMARY {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dde1efad90c922a0198cae64c642be95605c23cf1e53e3c35491817bf6c48b/data':
    'community.wave.seqera.io/library/pandas:3.1.0--f6c858a02e385828' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    path "*.csv", emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python ${projectDir}/modules/local/tRNA_summary/templates/gettRNAxlsites_v3.py -i ${bedfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS

    """
}
