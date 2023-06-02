/*
 * Check input samplesheet, get read channels and merge replicates if necessary.
 * Subworkflow goal is to output a valid, prepared set of fastq files with metadata
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check/main'
include { CAT_FASTQ         } from '../../modules/nf-core/cat/fastq/main'

workflow PARSE_FASTQ_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Check the samplesheet for errors
    */
    SAMPLESHEET_CHECK (
        samplesheet
    )

    /*
    * MODULE: Parse samplesheet into meta and fastq files
    */
    ch_reads = SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )
        .map { parse_meta(it) }

    /*
    * CHANNEL: Split out files which need merging
    */
    ch_fastq = ch_reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split("_")[0..-2].join("_")
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ (
        ch_fastq.multiple
    )
    ch_versions  = ch_versions.mix(CAT_FASTQ.out.versions)
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    emit:
    fastq    = ch_cat_fastq // channel: [ val(meta), [ reads ] ]
    versions = ch_versions  // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def parse_meta(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.id
    meta.group      = row.group
    meta.replicate  = row.replicate.toInteger()
    meta.single_end = row.single_end.toBoolean()

    // Get the rest of the values
    def keys = row.keySet()
    for (key in keys) {
        if (key in ['id', 'group', 'replicate', 'single_end', 'fastq_1', 'fastq_2']) continue
        meta[key] = row[key]
    }

    // Check fastq files exist
    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}
