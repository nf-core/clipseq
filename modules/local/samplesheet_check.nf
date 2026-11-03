process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container 'quay.io/biocontainers/pandas:1.5.2'

    input:
    path samplesheet
    val source

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/clipseq/bin/
    if (source == 'fastq') {
        """
        samplesheet_fastq_check.py --samplesheet $samplesheet --output samplesheet.valid.csv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    } else if (source == 'dedupe_bam') {
        """
        samplesheet_dedupe_bam_check.py --samplesheet $samplesheet --output samplesheet.valid.csv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}

