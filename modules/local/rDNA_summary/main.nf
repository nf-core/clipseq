process rDNA_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/4522cb2ea15e8a3a47ec2256b735bc6489807b9a444d01324ab722c652210f04/data':
    'community.wave.seqera.io/library/r-ggridges_r-tidyverse:5f6979fe527ac8bd' }"

    input:
    path(norm_bedgraph)

    output:
    path("*.png")                                    ,emit: png
    path  "versions.yml"                            ,emit: versions

    script:
    def args = task.ext.args ?: ''
    """

    Rscript $projectDir/modules/local/rDNA_summary/templates/rDNA_summary.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        r-ggridges: \$(Rscript -e "library(ggridges); cat(as.character(packageVersion('ggridges')))")
    END_VERSIONS
    """
}
