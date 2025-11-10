process bsfQC {
    tag "$meta.id"
    label 'process_single'

    container "${'docker.io/melinak/bindingsitefinder:1.1'}"

    input:
        tuple val(meta), path(binding_sites_rds)

    output:
        tuple val(meta), path("*BindingSiteFinderQC.html"), emit: bindingSiteFinderQC
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript BindingSiteFinderQC.R \\
        $binding_sites_rds
		
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