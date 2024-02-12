process LINUX_COMMAND {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta) , path(input)
    path input2
    val copy_input

    output:
    tuple val(meta), path("*.cmd.*"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def ext    = task.ext.ext ?: 'txt'
    def cmd1   = task.ext.cmd1 ?: 'echo "NO-ARGS"'
    def cmd2   = task.ext.cmd2 ? "CMD2=`cat $input2 | ${task.ext.cmd2}`" : ''
    if(copy_input) {
        cmd2 = task.ext.cmd2 ? "CMD2=`cat $input | ${task.ext.cmd2}`" : ''
    }

    """
    $cmd2
    cat $input | $cmd1 > ${prefix}.cmd.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linux: NOVERSION
    END_VERSIONS
    """
}
