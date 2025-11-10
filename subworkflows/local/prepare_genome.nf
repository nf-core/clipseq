// Import the log object
import nextflow.util.LoggerHelper

//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                                                 } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_NCRNA_FASTA                                           } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                                                   } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BT                                                      } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR                                                    } from '../../modules/nf-core/untar/main'
include { BOWTIE_BUILD                                                           } from '../../modules/nf-core/bowtie/build/main'
include { STAR_GENOMEGENERATE                                                    } from '../../modules/nf-core/star/genomegenerate/main'
include { SAMTOOLS_FAIDX as GENOME_INDEX                                         } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as NCRNA_INDEX                                          } from '../../modules/nf-core/samtools/faidx/main'
include { LINUX_COMMAND as REMOVE_GTF_BRACKETS                                   } from '../../modules/local/linux_command'
include { CUSTOM_GETCHROMSIZES as GENOME_CHROM_SIZE                              } from '../../modules/nf-core/custom/getchromsizes/main'
include { CUSTOM_GETCHROMSIZES as NCRNA_CHROM_SIZE                               } from '../../modules/nf-core/custom/getchromsizes/main'
include { FILTER_GTF_BY_TRANSCRIPT                                               } from '../../modules/local/filter_gtf_by_transcript/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_GTF                                   } from '../../modules/nf-core/icountmini/segment/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_FILTGTF                               } from '../../modules/nf-core/icountmini/segment/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_REGIONS             } from '../../modules/local/resolve_unannotated/main'

workflow PREPARE_GENOME {
    take:
    fasta                          // file: .fasta
    fasta_fai                      // file: .fai
    ncrna_fasta                    // file: .fasta
    ncrna_fasta_fai                // file: .fai
    gtf                            // file: .gtf
    genome_index                   // folder: index
    ncrna_genome_index             // folder: index
    genome_chrom_sizes             // file: .txt
    ncrna_chrom_sizes              // file: .txt
    representative_transcript      // file: .txt
    representative_transcript_fai  // file: .fai
    representative_transcript_gtf  // file: .gtf
    filtered_gtf                   // file: .gtf
    seg_gtf                        // file: .gtf
    regions_gtf                    // file: .gtf
    regions_filt_gtf               // file: .gtf
    regions_resolved_gtf           // file: .gtf
    skip_filter_gtf                // value: boolean
    skip_transcriptome             // value: boolean

    main:

    // Init
    ch_versions = Channel.empty()
    //
    // MODULE: Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.toString().endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( [ [id:fasta.baseName], fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.of([ [id:fasta.baseName], fasta ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fasta]
    //ch_fasta | view

    //
    // MODULE: Uncompress genome ncrna_fasta file if required
    //
    ch_ncrna_fasta = Channel.empty()
    if (ncrna_fasta.toString().endsWith(".gz")) {
        ch_ncrna_fasta = GUNZIP_NCRNA_FASTA ( [ [id:ncrna_fasta.baseName], ncrna_fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_NCRNA_FASTA.out.versions)
    } else {
        ch_ncrna_fasta = Channel.of([ [id:ncrna_fasta.baseName], ncrna_fasta ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fasta]
    //ch_ncrna_fasta | view

    //
    // MODULE: Uncompress genome gtf file if required
    //
    ch_gtf = Channel.empty()
    if (gtf.toString().endsWith(".gz")) {
        ch_gtf      = GUNZIP_GTF ( [ [id:gtf.baseName], gtf ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.of([ [id:gtf.baseName], gtf ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ch_gtf | view

    //
    // MODULES: Uncompress STAR index or generate if required
    //
    ch_star_index = Channel.empty()
    if (genome_index) {
        if (genome_index.toString().endsWith(".tar.gz")) {
            ch_star_index = UNTAR_STAR ( [ [:], genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_STAR.out.versions)
        } else {
            ch_star_index = Channel.of([ [:] , genome_index ])
        }
    }
    else {
        ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // MODULES: Uncompress Bowtie index or generate if required
    //
    ch_bt_index = Channel.empty()
    if (ncrna_genome_index) {
        if (ncrna_genome_index.toString().endsWith(".tar.gz")) {
            ch_bt_index = UNTAR_BT ( [ [:], ncrna_genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BT.out.versions)
        } else {
            ch_bt_index = Channel.of([ [:] , ncrna_genome_index ])
        }
    }
    else {
        ch_bt_index = BOWTIE_BUILD ( ch_ncrna_fasta ).index
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    }

    //
    // MODULE: Create fasta fai if required
    //
    ch_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_fasta_fai = Channel.of([ [id:fasta_fai.baseName], fasta_fai ])
    } else {
        GENOME_INDEX (
            ch_fasta,
            [[],[]],
	    false
        )
        ch_fasta_fai = GENOME_INDEX.out.fai
        ch_versions = ch_versions.mix(GENOME_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Create fasta fai if required for ncrna genome
    //
    ch_ncrna_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_ncrna_fasta_fai = Channel.of([ [id:ncrna_fasta_fai.baseName], fasta_fai ])
    } else {
        NCRNA_INDEX (
            ch_ncrna_fasta,
            [[],[]],
            false
        )
        ch_ncrna_fasta_fai = NCRNA_INDEX.out.fai
        ch_versions = ch_versions.mix(NCRNA_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Calc genome chrom sizes
    //
    ch_genome_chrom_sizes = Channel.empty()
    if(genome_chrom_sizes) {
        ch_genome_chrom_sizes = Channel.of([ [id:genome_chrom_sizes.baseName], genome_chrom_sizes ])
    } else {
        ch_genome_chrom_sizes = GENOME_CHROM_SIZE ( ch_fasta ).sizes
        ch_versions  = ch_versions.mix(GENOME_CHROM_SIZE.out.versions)
    }

    //
    // MODULE: Calc ncrna chrom sizes
    //
    ch_ncrna_chrom_sizes = Channel.empty()
    if(ncrna_chrom_sizes) {
        ch_ncrna_chrom_sizes = Channel.of([ [id:ncrna_chrom_sizes.baseName], ncrna_chrom_sizes ])
    } else {
        ch_ncrna_chrom_sizes = NCRNA_CHROM_SIZE ( ch_ncrna_fasta ).sizes
        ch_versions  = ch_versions.mix(NCRNA_CHROM_SIZE.out.versions)
    }

    //
    // MODULE: Remove brackets from in gene names from GTF as causes UMICollapse to fail.
    //
    REMOVE_GTF_BRACKETS (
        ch_gtf,
        [],
        false
    )
    ch_gtf = REMOVE_GTF_BRACKETS.out.file
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_gtf | view

    //
    // MODULE: Filter GTF by transcripts, validate GTF and transcripts (if provided), and generate auxiliary files (FAI, GTF) for provided or auto-selected transcripts
    //

    // Channels for representative_transcript
    ch_representative_transcript = representative_transcript ?
        Channel.of([ [id: representative_transcript.baseName], representative_transcript ]) :
        Channel.of([[], []])

    ch_representative_transcript_fai = representative_transcript_fai ?
        Channel.of([[id: representative_transcript_fai.baseName], representative_transcript_fai]) :
        Channel.empty()

    ch_representative_transcript_gtf = representative_transcript_gtf ?
        Channel.of([[id: representative_transcript_gtf.baseName], representative_transcript_gtf]) :
        Channel.empty()

    ch_filt_gtf = filtered_gtf ?
        Channel.of([[id: filtered_gtf.baseName], filtered_gtf]) :
        Channel.empty()

    // Run FILTER_GTF_BY_TRANSCRIPT if skip_filter_gtf disabled and filtered GTF missing,
    // OR any required transcript file (.txt, .fai, .gtf) is missing and skip_transcriptome false,
    // OR to validate the GTF and transcripts provided when representative_transcript provided and skip_transcriptome false, even when skip_filter_gtf is true.
    if (
        (!skip_filter_gtf && !filtered_gtf) ||
        ((!representative_transcript || !representative_transcript_fai || !representative_transcript_gtf) && !skip_transcriptome) ||
        ((representative_transcript && !skip_transcriptome && skip_filter_gtf))
    ) {
        if (representative_transcript && !skip_transcriptome && skip_filter_gtf) {
            log.info "INFO: You provided a representative_transcript file and enabled transcriptome analysis but set skip_filter_gtf=true. The FILTER_GTF_BY_TRANSCRIPT process will still run, only to validate your transcripts against the GTF."
        }
        FILTER_GTF_BY_TRANSCRIPT(ch_gtf, ch_representative_transcript, skip_filter_gtf)

        ch_representative_transcript     = FILTER_GTF_BY_TRANSCRIPT.out.representative_transcript
        ch_representative_transcript_fai = FILTER_GTF_BY_TRANSCRIPT.out.representative_transcript_fai
        ch_representative_transcript_gtf = FILTER_GTF_BY_TRANSCRIPT.out.representative_transcript_gtf
        if (!skip_filter_gtf && !filtered_gtf) {
            ch_filt_gtf                  = FILTER_GTF_BY_TRANSCRIPT.out.gtf
        }
        ch_versions                      = ch_versions.mix(FILTER_GTF_BY_TRANSCRIPT.out.versions)
    }

    //
    // MODULE: Segment GTF file using icount
    //
    ch_seg_gtf = seg_gtf ?
        Channel.of([[id: seg_gtf.baseName], seg_gtf]) :
        Channel.empty()

    ch_regions_gtf = regions_gtf ?
        Channel.of([[id: regions_gtf.baseName], regions_gtf]) :
        Channel.empty()

    if (!seg_gtf || !regions_gtf) {
        ICOUNT_SEG_GTF (
            ch_gtf,
            ch_fasta_fai.map{ it[1] }
        )
        ch_seg_gtf     = ICOUNT_SEG_GTF.out.gtf
        ch_regions_gtf = ICOUNT_SEG_GTF.out.regions
        ch_versions    = ch_versions.mix(ICOUNT_SEG_GTF.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ICOUNT_SEG_GTF.out.gtf | view

    //
    // MODULE: Segment the filtered GTF file using icount
    //
    ch_regions_filt_gtf = regions_filt_gtf ?
        Channel.of([[id: regions_filt_gtf.baseName], regions_filt_gtf]) :
        Channel.empty()

    if (!skip_filter_gtf && !regions_filt_gtf) {
        ICOUNT_SEG_FILTGTF (
            ch_filt_gtf,
            ch_fasta_fai.map{ it[1] }
        )

        ch_regions_filt_gtf = ICOUNT_SEG_FILTGTF.out.regions
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ICOUNT_SEG_FILTGTF.out.gtf | view


    //
    // MODULE: Resolve the GTF regions that iCount did not annotate REGIONS FILE
    //
    ch_regions_resolved_gtf = regions_resolved_gtf ?
        Channel.of([[id: regions_resolved_gtf.baseName], regions_resolved_gtf]) :
        Channel.empty()

    if (!skip_filter_gtf && !regions_resolved_gtf) {
        RESOLVE_UNANNOTATED_REGIONS (
            ch_regions_gtf,
            ch_regions_filt_gtf,
            ch_fasta_fai
        )
        ch_regions_resolved_gtf = RESOLVE_UNANNOTATED_REGIONS.out.gtf
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED_REGIONS.out.gtf | view



    emit:
    fasta                             = ch_fasta                             // channel: [ val(meta), [ fasta ] ]
    fasta_fai                         = ch_fasta_fai                         // channel: [ val(meta), [ fai ] ]
    ncrna_fasta                       = ch_ncrna_fasta                       // channel: [ val(meta), [ fasta ] ]
    ncrna_fasta_fai                   = ch_ncrna_fasta_fai                   // channel: [ val(meta), [ fai ] ]
    genome_index                      = ch_star_index                        // channel: [ val(meta), [ star_index ] ]
    ncrna_index                       = ch_bt_index                          // channel: [ val(meta), [ bt2_index ] ]
    chrom_sizes                       = ch_genome_chrom_sizes                // channel: [ val(meta), [ txt ] ]
    ncrna_chrom_sizes                 = ch_ncrna_chrom_sizes                 // channel: [ val(meta), [ txt ] ]
    gtf                               = ch_gtf                               // channel: [ val(meta), [ gtf ] ]
    representative_transcript         = ch_representative_transcript         // channel: [ val(meta), [ txt ] ] or Channel.empty()
    representative_transcript_fai     = ch_representative_transcript_fai     // channel: [ val(meta), [ fai ] ] or Channel.empty()
    representative_transcript_gtf     = ch_representative_transcript_gtf     // channel: [ val(meta), [ fai ] ] or Channel.empty()
    filtered_gtf                      = ch_filt_gtf                          // channel: [ val(meta), [ gtf ] ] or Channel.empty()
    seg_gtf                           = ch_seg_gtf                           // channel: [ val(meta), [ gtf ] ]
    regions_gtf                       = ch_regions_gtf                       // channel: [ val(meta), [ gtf ] ]
    regions_filt_gtf                  = ch_regions_filt_gtf                  // channel: [ val(meta), [ gtf ] ] or Channel.empty()
    regions_resolved_gtf              = ch_regions_resolved_gtf              // channel: [ val(meta), [ gtf ] ] or Channel.empty()
    versions                          = ch_versions                          // channel: [ versions.yml ]
}
