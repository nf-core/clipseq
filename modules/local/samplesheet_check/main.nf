process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.10.4"
    container "quay.io/biocontainers/python:3.10.4"

    input:
    path samplesheet
    val source

    output:
    path '*.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // bundled in bin folder
    process_name = task.process
    output       = task.ext.output ?: 'samplesheet.valid.csv'
    switch (source) {
        case 'fastq':
            """
            check_samplesheet_fastq_check.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
        case 'bam':
            """
            check_samplesheet_bam_check.py $samplesheet samplesheet.valid.csv
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
            END_VERSIONS
            """
            break;
    }
}
