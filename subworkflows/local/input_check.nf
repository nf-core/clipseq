include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    source      // value: params.source

    main:

    ch_versions = Channel.empty() // initialise for this workflow

    switch(source) {
        case 'fastq': // fastq ready for alignment
            SAMPLESHEET_CHECK ( samplesheet, source )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .groupTuple(by: [0])
            .branch {
                meta, fastq ->
                    single  : fastq.size() == 1
                        return [ meta, fastq.flatten() ]
                    multiple: fastq.size() > 1
                        return [ meta, fastq.flatten() ]
            }
             .set { ch_fastq }
             CAT_FASTQ (
                 ch_fastq.multiple
             )
             ch_versions  = ch_versions.mix(CAT_FASTQ.out.versions)
             reads = CAT_FASTQ.out.reads.mix(ch_fastq.single)
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
def create_dedupe_bam_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id            = row.sample_name
    meta.group         = row.group_name
    meta.control       = row.input_name

    // Check fastq files exist
    def array = []
    if (!file(row.dedupe_bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 bam file does not exist!\n${row.fastq}"
    }
    array = [ meta, [ file(row.dedupe_bam) ] ]
    return array

    // Check if grouping needs to happen
}
