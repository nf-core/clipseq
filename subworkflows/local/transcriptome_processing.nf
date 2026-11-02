include { SAMTOOLS_SIMPLE_VIEW as FILTER_TRANSCRIPTS                          } from '../../modules/local/samtools_simple_view'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as TRANS_DEDUP                       } from './bam_dedup_samtools_umicollapse'
include { GET_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL             } from '../../modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL_HASGROUP    } from '../../modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS_GROUP_HASGROUP         } from '../../modules/local/get_crosslinks'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILT_TRANSCRIPT                      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILT_TRANSCRIPT                    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_TRANSCRIPTOME                      } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_TRANSCRIPTOME_GROUP_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { CLIPPY as CLIPPY_TRANSCRIPTOME                                      } from "../../modules/nf-core/clippy/main"
include { PARACLU as PARACLU_TRANSCRIPTOME                                    } from "../../modules/nf-core/paraclu/main"



workflow TRANSCRIPTOME_PROCESSING {
    take:
    ch_transcript_bam          // channel: [ val(meta), [ bam ] ]
    ch_transcript_bai          // channel: [ val(meta), [ bai ] ]
    ch_representative_transcript      // channel: []
    ch_representative_transcript_gtf  // channel: []
    ch_representative_transcript_fai  // channel: []
    callers
    ch_paraclu_mincluster

    main:
    ch_versions = Channel.empty()
    if(params.source == "fastq" & params.run_filtering) {
        //
        // CHANNEL: Combine bam and bai files on id
        //
        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }
        //ch_transcript_bam_bai | view

        //
        // MODULE: Filter transcriptome bam on longest transcripts
        //
        FILTER_TRANSCRIPTS (
            ch_transcript_bam_bai,
            [],
            ch_representative_transcript
        )
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS.out.versions)

        //
        // Sort, index filtered bam
        //
        SAMTOOLS_SORT_FILT_TRANSCRIPT ( FILTER_TRANSCRIPTS.out.bam , [[],[]] )
        ch_versions       = ch_versions.mix(SAMTOOLS_SORT_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bam = SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam

        SAMTOOLS_INDEX_FILT_TRANSCRIPT ( SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam )
        ch_versions       = ch_versions.mix(SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bai = SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.bai
    }

    // UMI DEDUPLICATION
    ch_trans_umi_log  = Channel.empty()
    if(params.source == "fastq" & params.run_dedup) {
        ch_transcript_bam_bai = ch_transcript_bam
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }

        //
        // SUBWORKFLOW: Run umi deduplication on transcript-level alignments
        //
        TRANS_DEDUP (
            ch_transcript_bam_bai
        )
        ch_versions        = ch_versions.mix(TRANS_DEDUP.out.versions)
        ch_transcript_bam  = TRANS_DEDUP.out.bam
        ch_transcript_bai  = TRANS_DEDUP.out.bai
        ch_trans_umi_log   = TRANS_DEDUP.out.umi_log
    }

    ch_trans_crosslink_bed            = Channel.empty()
    ch_trans_crosslink_coverage       = Channel.empty()
    ch_trans_crosslink_coverage_norm  = Channel.empty()

    //
    // GROUPING: At this point, if groups have been specified, then we need to merge corresponding BAM files
    //
    // ch_genome_bam.view { item -> println("Pre-branch: $item") }

    ch_transcript_bam.branch {
        hasGroup: it[0].group  // Branch condition for samples with a group
        noGroup: it[0].group == ''  // Branch condition for samples without a group
    }.set { ch_bam_branches }  // Capture branching result into ch_branches

    ch_transcript_bai.branch {
        hasGroup: it[0].group  // Branch condition for samples with a group
        noGroup: it[0].group == ''  // Branch condition for samples without a group
    }.set { ch_bai_branches }  // Capture branching result into ch_branches

    // After branching
    //ch_branches.hasGroup.view { item -> "Has group: $item" }
    //ch_branches.noGroup.view { item -> "No group: $item" }

    ch_bam_branches.hasGroup
        .map { item ->
            def meta = item[0]
            def bam = item[1]
            return [meta.group, meta, bam]
        }
        .groupTuple(by: 0)
        .map { tuple ->
            def group = tuple[0]
            def items = tuple[1]
            def bam = tuple[2]

            def newMeta = [:]
            newMeta.id = group
            newMeta.group = group
            newMeta.control = items[0].control
            newMeta.single_end = true

            return [newMeta, bam]
        }
        //.view { "Grouped and remapped: $it" }
        .set { ch_grouped_transcript_bam }

    SAMTOOLS_MERGE_TRANSCRIPTOME (
        ch_grouped_transcript_bam,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_TRANSCRIPTOME.out.versions)
    ch_grouped_transcript_bam = SAMTOOLS_MERGE_TRANSCRIPTOME.out.bam

    //ch_grouped_transcript_bam.view { item -> "Grouped transcript BAM: $item" }

    // merge with the noGroup branch
    ch_grouped_transcript_bam.concat(ch_bam_branches.noGroup).set { ch_transcript_peakcalling_bam }
    //ch_genome_peakcalling_bam.view { item -> "Grouped genome BAM and non-grouped BAMs: $item" }

    SAMTOOLS_TRANSCRIPTOME_GROUP_INDEX (
        ch_transcript_peakcalling_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_TRANSCRIPTOME_GROUP_INDEX.out.versions)
    ch_transcript_peakcalling_bai = SAMTOOLS_TRANSCRIPTOME_GROUP_INDEX.out.bai

    if(params.run_crosslinking) {
        //
        // SUBWORKFLOW: Run crosslink calculation for transcripts
        //
        CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL (
            ch_bam_branches.noGroup.join(ch_bai_branches.noGroup),
            ch_representative_transcript_fai.map{ [[id:it.baseName], it] }
        )
        ch_versions                       = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL.out.versions)
        ch_trans_crosslink_INDIVIDUAL_bed = CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL.out.bed

        CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL_HASGROUP (
            ch_bam_branches.hasGroup.join(ch_bai_branches.hasGroup),
            ch_representative_transcript_fai.map{ [[id:it.baseName], it] }
        )
        ch_versions                                   = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL_HASGROUP.out.versions)
        ch_trans_crosslink_INDIVIDUAL_HASGROUP_bed    = CALC_TRANSCRIPT_CROSSLINKS_INDIVIDUAL_HASGROUP.out.bed

        CALC_TRANSCRIPT_CROSSLINKS_GROUP_HASGROUP (
            ch_grouped_transcript_bam.join(ch_transcript_peakcalling_bai),
            ch_representative_transcript_fai.map{ [[id:it.baseName], it] }
        )
        ch_versions                            = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS_GROUP_HASGROUP.out.versions)
        ch_trans_crosslink_GROUP_HASGROUP_bed  = CALC_TRANSCRIPT_CROSSLINKS_GROUP_HASGROUP.out.bed

        //
        // Combine crosslinking results for moving forwards
        //
        ch_trans_crosslink_bed = ch_trans_crosslink_INDIVIDUAL_bed.mix(ch_trans_crosslink_GROUP_HASGROUP_bed)

    }


    ch_clippy_transcriptome_peaks   = Channel.empty()
    ch_paraclu_transcriptome_peaks  = Channel.empty()
    if(params.run_peakcalling) {
        if('clippy' in callers) {
            CLIPPY_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_representative_transcript_gtf,
                ch_representative_transcript_fai
            )

            ch_clippy_transcriptome_peaks    = CLIPPY_TRANSCRIPTOME.out.peaks
            ch_versions                      = ch_versions.mix(CLIPPY_TRANSCRIPTOME.out.versions)
        }

        if('paraclu' in callers) {
            PARACLU_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_paraclu_mincluster
            )
            ch_paraclu_transcriptome_peaks          = PARACLU_TRANSCRIPTOME.out.bed
            ch_versions                             = ch_versions.mix(PARACLU_TRANSCRIPTOME.out.versions)
        }
    }


    emit:
    transcript_dedupe_bam   = ch_transcript_peakcalling_bam      // channel: [ val(meta), [ bam ] ]
    transcript_dedupe_bai   = ch_transcript_peakcalling_bai      // channel: [ val(meta), [ bai ] ]

    crosslink_bed           = ch_trans_crosslink_bed             // channel: [ val(meta), [ bed ] ]

    clippy_peaks            = ch_clippy_transcriptome_peaks      // channel: [ val(meta), [ bed ] ]
    paraclu_peaks           = ch_paraclu_transcriptome_peaks     // channel: [ val(meta), [ bed ] ]

    versions                = ch_versions                        // channel: [ versions.yml ]
}
