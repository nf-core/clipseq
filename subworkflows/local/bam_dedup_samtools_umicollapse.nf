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
    fasta   // channel: [ val(meta), [ fasta/fa ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: UMI-tools collapse
    //
    UMICOLLAPSE (
        bam_bai 
    )
    ch_versions = ch_versions.mix(UMICOLLAPSE.out.versions)

    //
    // MODULE: Index BAM file
    //
    SAMTOOLS_INDEX (
        UMICOLLAPSE.out.bam 
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

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

    //
    // SUBWORKFLOW: Re-calc stats
    //
    BAM_STATS_SAMTOOLS ( ch_bam_bai, fasta )

    emit:
    bam      = UMICOLLAPSE.out.bam             // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    umi_log  = UMICOLLAPSE.out.log             // channel: [ val(meta), [ log ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
