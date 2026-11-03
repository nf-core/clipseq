//
// UMIcollapse, index BAM file
//

//
// MODULES
//

include { UMICOLLAPSE    } from '../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

//
// SUBWORKFLOWS
//

include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow BAM_DEDUP_SAMTOOLS_UMICOLLAPSE {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai/csi ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: UMI-tools collapse
    //
    UMICOLLAPSE (
        bam_bai,
        'bam'
    )

    //
    // MODULE: Index BAM file
    //
    SAMTOOLS_INDEX (
        UMICOLLAPSE.out.bam
    )

    //
    // CHANNEL: Merge reads and index
    //
    ch_bam_bai = UMICOLLAPSE.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }


    emit:
    bam      = UMICOLLAPSE.out.bam             // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    umi_log  = UMICOLLAPSE.out.log             // channel: [ val(meta), [ log ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
