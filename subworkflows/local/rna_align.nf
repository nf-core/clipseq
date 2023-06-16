//
// Align to smrna and primary genomes before indexing and sorting
//

//
// MODULES
//
include { BOWTIE_ALIGN } from '../../modules/nf-core/bowtie/align/main.nf'
include { STAR_ALIGN   } from '../../modules/nf-core/star/align/main.nf'

//
// SUBWORKFLOWS
//
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_SMRNA             } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_TARGET_GENOME     } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_TARGET_TRANSCRIPT } from '../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow RNA_ALIGN {
    take:
    fastq      // channel: [ val(meta), [ fastq ] ]
    bt2_index  // channel: [ val(meta), [ index ] ]
    star_index // channel: [ val(meta), [ index ] ]
    gtf        // channel: [ val(meta), [ gtf ] ]
    fasta      // channel: [ val(meta), [ fasta/fa ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Align reads to smrna genome
    //
    BOWTIE_ALIGN (
        fastq,
        bt2_index.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS_SMRNA ( BOWTIE_ALIGN.out.bam, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_SMRNA.out.versions)

    //
    // MODULE: Align reads that did not align to the smrna genome to the primary genome
    //
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        star_index.map { it[1] },
        gtf.map { it[1] },
        false,
        '',
        ''
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS_TARGET_GENOME ( STAR_ALIGN.out.bam_sorted, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_TARGET_GENOME.out.versions)

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS_TARGET_TRANSCRIPT ( STAR_ALIGN.out.bam_transcript, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_TARGET_TRANSCRIPT.out.versions)

    // emit:
    // bt_bam              = SAMTOOLS_SORT_SMRNA.out.bam                            // channel: [ val(meta), [ bam ] ]
    // bt_bai              = SAMTOOLS_INDEX_SMRNA.out.bai                           // channel: [ val(meta), [ bam ] ]
    // bt_log              = BOWTIE_ALIGN.out.log                            // channel: [ val(meta), [ txt ] ]
    // star_bam            = STAR_ALIGN.out.bam_sorted                       // channel: [ val(meta), [ bam ] ]
    // star_log            = STAR_ALIGN.out.log                              // channel: [ val(meta), [ txt ] ]
    // star_log_final      = STAR_ALIGN.out.log_final                        // channel: [ val(meta), [ txt ] ]
    // genome_bam          = STAR_ALIGN.out.bam_sorted                       // channel: [ val(meta), [ bam ] ]
    // genome_bai          = SAMTOOLS_INDEX_GENOME.out.bai                          // channel: [ val(meta), [ bai ] ]
    // transcript_bam = SAMTOOLS_SORT_TRANSCRIPT.out.bam      // channel: [ val(meta), [ bam ] ]
    // transcript_bai      = SAMTOOLS_INDEX_TRANSCRIPT.out.bai      // channel: [ val(meta), [ bai ] ]
    // versions            = ch_versions                                     // channel: [ versions.yml ]
}
