//
// Calculate clip crosslinks using an input BAM file and genome index file.
// Crosslinks are outputed as a BED file and additional coverage and normalised coverage
// tracks are calculated in BEDGRAPH format
//

//
// MODULES
//
include { BEDTOOLS_BAMTOBED                            } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_SHIFT                               } from '../../modules//nf-core/bedtools/shift/main'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS } from '../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG } from '../../modules/nf-core/bedtools/genomecov/main'
include { LINUX_COMMAND as SELECT_BED_POS              } from '../../modules/local/linux_command'
include { LINUX_COMMAND as SELECT_BED_NEG              } from '../../modules/local/linux_command'
include { LINUX_COMMAND as MERGE_AND_SORT              } from '../../modules/local/linux_command'
include { LINUX_COMMAND as CROSSLINK_COVERAGE          } from '../../modules/local/linux_command'
include { LINUX_COMMAND as CROSSLINK_NORMCOVERAGE      } from '../../modules/local/linux_command'

workflow CALC_CROSSLINKS {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    fai // channel: [ fai ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Convert input BAM file to BED format
    //
    BEDTOOLS_BAMTOBED (
        bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    //
    // MODULE: Shift BED file according to parameters suppied in config (default is -s 0)
    //
    BEDTOOLS_SHIFT (
        BEDTOOLS_BAMTOBED.out.bed,
        fai
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SHIFT.out.versions)

    //
    // MODULE: Report depth at each position on the pos strand
    //
    BEDTOOLS_GENOMECOV_POS (
        BEDTOOLS_SHIFT.out.bed.map{ [ it[0], it[1], 1 ] },
        fai.map{ it[1] },
        'pos.bed'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_POS.out.versions)

    //
    // MODULE: Report depth at each position on the neg strand
    //
    BEDTOOLS_GENOMECOV_NEG (
        BEDTOOLS_SHIFT.out.bed.map{ [ it[0], it[1], 1 ] },
        fai.map{ it[1] },
        'neg.bed'
    )

    //
    // MODULE: Select columns in BED file using AWK
    //
    SELECT_BED_POS (
        BEDTOOLS_GENOMECOV_POS.out.genomecov,
        [],
        false
    )
    SELECT_BED_NEG (
        BEDTOOLS_GENOMECOV_NEG.out.genomecov,
        [],
        false
    )

    //
    // CHANNEL: Join POS/NEG files into one channel so they can be merged in the next module
    //
    ch_merge_and_sort_input = SELECT_BED_POS.out.file
        .map{ [ it[0].id, it[0], it[1] ] }
        .join( SELECT_BED_NEG.out.file.map{ [ it[0].id, it[0], it[1] ] } )
        .map { [ it[1], [ it[2], it[4] ] ] }
    //EXAMPLE CHANNEL STRUCT: [ [id:test], [ BED(pos), BED(neg) ] ]
    //ch_merge_and_sort_input | view 

    //
    // MODULE: Select columns in BED file using AWK
    //
    MERGE_AND_SORT (
        ch_merge_and_sort_input,
        [],
        false
    )

    //
    // MODULE: Create coverage track using AWK
    //
    CROSSLINK_COVERAGE (
        MERGE_AND_SORT.out.file,
        [],
        false
    )

    //
    // MODULE: Create normalised coverage track using AWK
    //
    CROSSLINK_NORMCOVERAGE (
        MERGE_AND_SORT.out.file,
        [],
        true
    )

    emit:
    bed           = MERGE_AND_SORT.out.file         // channel: [ val(meta), [ bed ] ]
    coverage      = CROSSLINK_COVERAGE.out.file     // channel: [ val(meta), [ bedgraph ] ]
    coverage_norm = CROSSLINK_NORMCOVERAGE.out.file // channel: [ val(meta), [ bedgraph ] ]
    versions      = ch_versions                     // channel: [ versions.yml ]
}
