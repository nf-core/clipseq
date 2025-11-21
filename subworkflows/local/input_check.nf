include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    source      // value: params.source

    main:

    switch(source) {
        case 'fastq': // fastq ready for alignment
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { reads }
            break;
        case 'bam': // bam ready for alignment
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_bam_channel(it) }
            .set { reads }
            break;
        case 'dedupe_bam': // dedupe bam ready for grouped crosslink/peak analysis
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_dedupe_bam_channel(it) }
            .set { reads }
            break;
    }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

// Function to get list of [ meta, [ fastq ] ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample_name
    meta.group      = row.group_name
    meta.control    = row.input_name
    meta.single_end = true

    // Check fastq files exist
    def array = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq}"
    }
    array = [ meta, [ file(row.fastq) ] ]
    return array

}


// Function to get list of [ meta, bam ]
def create_bam_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id            = row.sample_name
    meta.group         = row.group_name
    meta.control       = row.input_name

    // Check bam files exist
    def array = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
    }
    array = [ meta, [ file(row.bam) ] ]
    return array
}

// Function to get list of [ meta, bam ]
def create_dedupe_bam_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id            = row.sample_name
    meta.group         = row.group_name
    meta.control       = row.input_name

    // Check dedupe_bam files exist
    def array = []
    if (!file(row.dedupe_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Dedupe BAM file does not exist!\n${row.dedupe_bam}"
    }
    array = [ meta, [ file(row.dedupe_bam) ] ]
    return array

    // Check if grouping needs to happen
}
