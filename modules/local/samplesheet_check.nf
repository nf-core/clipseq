process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'
    errorStrategy 'terminate'
    debug true

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    val source

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/clipseq/bin/

    switch (source) {
        case 'fastq':
            """
            samplesheet_fastq_check.py --samplesheet $samplesheet --output samplesheet.valid.csv
            echo "Python script output:"
            cat .command.out
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
        case 'dedupe_bam':
            """
            samplesheet_dedupe_bam_check.py --samplesheet $samplesheet --output samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
    }
}

