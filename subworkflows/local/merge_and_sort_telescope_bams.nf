include { SAMTOOLS_VIEW as FILTER_UNIQUE_MAP_UPDATED                } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as FILTER_UNIQUE_MAP_OTHER                  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_MERGE as MERGE_TE_BAMS                           } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT as SORT_TE_BAMS                           } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as INDEX_TE_BAMS                           } from '../../modules/nf-core/samtools/index/main'

workflow MERGE_AND_SORT_TELESCOPE_BAMS {
    take:
    telescope_out_updated
    telescope_out_other

    main:
    // FILTER FOR UNIQUE MAPPERS ONLY FROM TELESCOPE OUTPUT

    FILTER_UNIQUE_MAP_UPDATED(
        telescope_out_updated.map { meta, bam -> [meta, bam, []] },
        [[],[]],
        []
    )

    FILTER_UNIQUE_MAP_OTHER(
        telescope_out_other.map { meta, bam -> [meta, bam, []] },
        [[],[]],
        []
        )

    // MERGE TELESCOPE BAMS TOGETHER FOR PEAK CALLING

    bams_to_merge = FILTER_UNIQUE_MAP_UPDATED.out.bam
        .join(FILTER_UNIQUE_MAP_OTHER.out.bam)
        .map { meta, bam1, bam2 -> [meta, [bam1, bam2]] }

    MERGE_TE_BAMS(
        bams_to_merge,
        [[],[]],
        [[],[]]
    )

    SORT_TE_BAMS(
        MERGE_TE_BAMS.out.bam,
        [[],[]]
    )

    INDEX_TE_BAMS(
        SORT_TE_BAMS.out.bam
    )

    emit:
    bam = SORT_TE_BAMS.out.bam
    bai = INDEX_TE_BAMS.out.bai
}