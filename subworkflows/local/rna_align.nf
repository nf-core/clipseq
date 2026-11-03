//
// Align to rDNA and primary genomes before indexing and sorting
//

//
// MODULES
//
include { BOWTIE_ALIGN                                      } from '../../modules/nf-core/bowtie/align/main.nf'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_K1                   } from '../../modules/nf-core/bowtie/align/main.nf'
include { STAR_ALIGN as STAR_ALIGN_WITH_TRANSCRIPTOME       } from '../../modules/nf-core/star/align/main.nf'
include { STAR_ALIGN as STAR_ALIGN_GENOME_ONLY              } from '../../modules/nf-core/star/align/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_COORD_TONLY      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_COORD_GONLY      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANS              } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANS            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_rDNA               } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_rDNA             } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_rDNA_K1            } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_rDNA_K1          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as TONLY_SAMTOOLS_VIEW_GENOME       } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as TONLY_SAMTOOLS_VIEW_TRANSCRIPT   } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as GONLY_SAMTOOLS_VIEW_GENOME       } from '../../modules/nf-core/samtools/view/main'

//
// SUBWORKFLOWS
//

workflow RNA_ALIGN {
    take:
    fastq               // channel: [ val(meta), [ fastq ] ]
    bt_index           // channel: [ val(meta), [ index ] ]
    star_index          // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    fasta               // channel: [ val(meta), [ fasta/fa ]
    skip_transcriptome  // boolean

    main:
    ch_versions = channel.empty()

    // SAMTOOLS_SORT expects [ meta, fasta, fai ]
    ch_sort_fasta = fasta.map { meta, fa -> [ meta, fa, [] ] }
    //
    // MODULE: Align reads to rDNA genome
    //
    BOWTIE_ALIGN (
        fastq,
        bt_index,
        true
    )

    //
    // SUBWORKFLOW: Sort, index BAM file
    //
    SAMTOOLS_SORT_rDNA( BOWTIE_ALIGN.out.bam, ch_sort_fasta, '' )
    SAMTOOLS_INDEX_rDNA( SAMTOOLS_SORT_rDNA.out.bam )



    /*
    * MODULE: Align reads to smrna genome, here allowing 100 multimappers but only reporting one alignment per multimapped read
    * so that we can accurately count it in the crosslink summary later
    */

    BOWTIE_ALIGN_K1 (
        fastq,
        bt_index,
        true
    )

    SAMTOOLS_SORT_rDNA_K1 ( BOWTIE_ALIGN_K1.out.bam, ch_sort_fasta, '' )

    SAMTOOLS_INDEX_rDNA_K1 ( SAMTOOLS_SORT_rDNA_K1.out.bam )

    //
    // MODULE: Align reads that did not align to the rDNA genome to the primary genome
    //
    if (skip_transcriptome) {
        STAR_ALIGN_GENOME_ONLY (
            BOWTIE_ALIGN.out.fastq,
            star_index,
            gtf,
            false
        )

        //
        // MODULE: Index the coord reads
        //
        SAMTOOLS_INDEX_COORD_GONLY( STAR_ALIGN_GENOME_ONLY.out.bam_sorted )

        //
        // CHANNEL: Merge reads and index
        //
        ch_coord_bam_bai = STAR_ALIGN_GENOME_ONLY.out.bam_sorted
            .join(SAMTOOLS_INDEX_COORD_GONLY.out.index, by: [0])
        /*
        * CHANNEL: Filter for uniquely mapping reads for downstream analysis; samtools view -b -q 5 -o output.bam alignments.bam
        */
        GONLY_SAMTOOLS_VIEW_GENOME (
            ch_coord_bam_bai,
            [[],[],[]],
            [[],[]],
            [[],[]],
            'bai'
        )

        ch_genome_log             = STAR_ALIGN_GENOME_ONLY.out.log_final
        ch_genome_log_final       = STAR_ALIGN_GENOME_ONLY.out.log_final
        ch_genome_multi_bam       = STAR_ALIGN_GENOME_ONLY.out.bam_sorted
        ch_genome_multi_bai       = SAMTOOLS_INDEX_COORD_GONLY.out.index
        ch_genome_unique_bam      = GONLY_SAMTOOLS_VIEW_GENOME.out.bam
        ch_genome_unique_bai      = GONLY_SAMTOOLS_VIEW_GENOME.out.bai
        ch_transcript_unique_bam  = []
        ch_transcript_unique_bai  = []
        ch_transcript_multi_bam   = []
        ch_transcript_multi_bai   = []
    } else {
        STAR_ALIGN_WITH_TRANSCRIPTOME (
            BOWTIE_ALIGN.out.fastq,
            star_index,
            gtf,
            false
        )

        //
        // MODULE: Index the coord reads
        //
        SAMTOOLS_INDEX_COORD_TONLY( STAR_ALIGN_WITH_TRANSCRIPTOME.out.bam_sorted )

        //
        // CHANNEL: Merge reads and index
        //
        ch_coord_bam_bai = STAR_ALIGN_WITH_TRANSCRIPTOME.out.bam_sorted
            .join(SAMTOOLS_INDEX_COORD_TONLY.out.index, by: [0])
        //
        // MODULE: Sort and index transcript BAM file
        //
        SAMTOOLS_SORT_TRANS( STAR_ALIGN_WITH_TRANSCRIPTOME.out.bam_transcript, ch_sort_fasta, '' )
        SAMTOOLS_INDEX_TRANS( SAMTOOLS_SORT_TRANS.out.bam )

        /*
        * CHANNEL: Join transcriptbam and bai files
        */
        ch_transcript_bam_bai = SAMTOOLS_SORT_TRANS.out.bam
            .join(SAMTOOLS_INDEX_TRANS.out.index, by: [0])
        /*
        * CHANNEL: Filter for uniquely mapping reads for downstream analysis; samtools view -b -q 5 -o output.bam alignments.bam
        */
        TONLY_SAMTOOLS_VIEW_GENOME (
            ch_coord_bam_bai,
            [[],[],[]],
            [[],[]],
            [[],[]],
            'bai'
        )

        TONLY_SAMTOOLS_VIEW_TRANSCRIPT (
            ch_transcript_bam_bai,
            [[],[],[]],
            [[],[]],
            [[],[]],
            'bai'
        )

        ch_genome_log             = STAR_ALIGN_WITH_TRANSCRIPTOME.out.log_final
        ch_genome_log_final       = STAR_ALIGN_WITH_TRANSCRIPTOME.out.log_final
        ch_genome_multi_bam       = STAR_ALIGN_WITH_TRANSCRIPTOME.out.bam_sorted
        ch_genome_multi_bai       = SAMTOOLS_INDEX_COORD_TONLY.out.index
        ch_genome_unique_bam      = TONLY_SAMTOOLS_VIEW_GENOME.out.bam
        ch_genome_unique_bai      = TONLY_SAMTOOLS_VIEW_GENOME.out.bai

        ch_transcript_unique_bam  = TONLY_SAMTOOLS_VIEW_TRANSCRIPT.out.bam
        ch_transcript_unique_bai  = TONLY_SAMTOOLS_VIEW_TRANSCRIPT.out.bai
        ch_transcript_multi_bam   = SAMTOOLS_SORT_TRANS.out.bam
        ch_transcript_multi_bai   = SAMTOOLS_INDEX_TRANS.out.index
    }


    emit:
    rDNA_bam         = SAMTOOLS_SORT_rDNA.out.bam       // channel: [ val(meta), [ bam ] ]
    rDNA_bai         = SAMTOOLS_INDEX_rDNA.out.index    // channel: [ val(meta), [ bai ] ]
    rDNA_log         = BOWTIE_ALIGN.out.log             // channel: [ val(meta), [ txt ] ]
    rDNA_k1_bam      = SAMTOOLS_SORT_rDNA_K1.out.bam                   // channel: [ val(meta), [ bam ] ]
    rDNA_k1_bai      = SAMTOOLS_INDEX_rDNA_K1.out.index                // channel: [ val(meta), [ bai ] ]

    genome_log             = ch_genome_log                    // channel: [ val(meta), [ txt ] ]
    genome_log_final       = ch_genome_log_final              // channel: [ val(meta), [ txt ] ]
    genome_unique_bam      = ch_genome_unique_bam                    // channel: [ val(meta), [ bam ] ]
    genome_unique_bai      = ch_genome_unique_bai                    // channel: [ val(meta), [ bai ] ]
    genome_multi_bam       = ch_genome_multi_bam                    // channel: [ val(meta), [ bam ] ]
    genome_multi_bai       = ch_genome_multi_bai                    // channel: [ val(meta), [ bai ] ]

    transcript_unique_bam  = ch_transcript_unique_bam                // channel: [ val(meta), [ bam ] ]
    transcript_unique_bai  = ch_transcript_unique_bai                // channel: [ val(meta), [ bai ] ]
    transcript_multi_bam   = ch_transcript_multi_bam                // channel: [ val(meta), [ bam ] ]
    transcript_multi_bai   = ch_transcript_multi_bai                // channel: [ val(meta), [ bai ] ]

    versions         = ch_versions                      // channel: [ versions.yml ]
}
