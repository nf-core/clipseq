process bsfQC {
    tag "$meta.id"
    label 'process_low'

        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://melinak/bindingsitefinder:latest':
        'melinak/bindingsitefinder:latest' }"

    input:
        tuple val(meta), path(binding_sites_rds)

    output:
        tuple val(meta), path("*BindingSiteFinderQC.html"), emit: BindingSiteFinderQC
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript clipseq/modules/local/BindingSiteFinder/bsfQC/BindingSiteFinderQC.R \\
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