process GET_INIT_ALIGNED_XLINKS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96dde1efad90c922a0198cae64c642be95605c23cf1e53e3c35491817bf6c48b/data':
    'community.wave.seqera.io/library/bedtools_pybedtools_pysam_matplotlib_pruned:79786472f5bd377f' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_init_xl_coord.bam")    , emit: bam
    path  "versions.yml"                            ,emit: versions

    script:
    template 'get_xlinks.py'
}
