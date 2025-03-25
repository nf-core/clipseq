process sortAnnotationForBindingSiteFinder {
    tag "$meta.id"
    label 'process_single'

    container "${'melinak/bindingsitefinder:1.1'}"

    input:
        tuple val(meta), path(gtf_file)

    output:
        tuple val(meta), path("*gns.rds"), emit: gns_rds
        tuple val(meta), path("*regions.rds"), emit: regions_rds
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript /home/mek24iv/nfcore-clipseq/devel_BindingSiteFinder/modules/local/sortAnnotationForBindingSiteFinder/sortAnnotationForBindingSiteFinder.R \\
        $gtf_file \\
        gns.rds \\
        regions.rds
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(Rscript -e "packageVersion('BindingSiteFinder')" |& sed '1!d ; s/[1]  //')
    END_VERSIONS
        """

    stub:
    def args = task.ext.a
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gns.rds
    touch regions.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(Rscript -e "packageVersion('BindingSiteFinder')" |& sed '1!d ; s/[1]  //')
    END_VERSIONS
    """
}