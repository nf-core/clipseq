process COMBINE_BINS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9b82e07dc06620144a486ca8e8b326c4594b4bc3563114235fb2ab377ee4093d/data':
        'community.wave.seqera.io/library/pandas_pyyaml:70aa423a6e4688f5' }"

    input:
    tuple val(meta) , path(bin_ncRNA)
    tuple val(meta2), path(bin_regions)
    path(premap_logs)

    output:
    tuple val(meta), path("*multimap_binning.csv"),  emit: csv
    path "colours_mqc_config.yaml", emit: config
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

     python ${projectDir}/modules/local/combine_bins/templates/combine_bins.py \
     -a ${bin_ncRNA} \
     -b ${bin_regions} \
     -premap ${premap_logs}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
