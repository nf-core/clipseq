include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    source      // value: params.source

    main:

    switch(source) {
        case 'fastq': // fastq ready for alignment
            println "samplesheet: $samplesheet, source: $source" // Print the inputs
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .view { item -> println "SAMPLESHEET_CHECK output: $item" } // Print the output
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { reads }
            println "Reads: $reads" // Print the reads
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
    println "Row: $row" // Print the row
    def meta = [:]
    meta.id         = row.sample_name
    meta.group      = row.group_name
    meta.control    = row.input_name
    meta.fastq      = row.fastq
    println "Meta: $meta" // Print the meta
    // Check fastq files exist
    def array = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq}"
    }
    array = [ meta, [ file(row.fastq) ] ]
    return array
}


// Function to get list of [ meta, bam ]
def create_dedupe_bam_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id            = row.sample_name
    meta.group         = row.group_name
    meta.control       = row.input_name
    meta.dedupe_bam    = row.dedupe_bam

    // Check fastq files exist
    def array = []
    if (!file(row.dedupe_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq}"
    }
    array = [ meta, [ file(row.dedupe_bam) ] ]
    return array
}
