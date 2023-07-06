//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                                                 } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SMRNA_FASTA                                           } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                                                   } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BT                                                      } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR                                                    } from '../../modules/nf-core/untar/main'
include { BOWTIE_BUILD                                                           } from '../../modules/nf-core/bowtie/build/main'
include { STAR_GENOMEGENERATE                                                    } from '../../modules/nf-core/star/genomegenerate/main'
include { SAMTOOLS_FAIDX as TARGET_INDEX                                         } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SMRNA_INDEX                                          } from '../../modules/nf-core/samtools/faidx/main'
include { LINUX_COMMAND as REMOVE_GTF_BRACKETS                                   } from '../../modules/local/linux_command'
include { CUSTOM_GETCHROMSIZES as TARGET_CHROM_SIZE                              } from '../../modules/nf-core/custom/getchromsizes/main'
include { CUSTOM_GETCHROMSIZES as SMRNA_CHROM_SIZE                               } from '../../modules/nf-core/custom/getchromsizes/main'
include { FIND_LONGEST_TRANSCRIPT                                                } from '../../modules/local/find_longest_transcript/main'
include { CLIPSEQ_FILTER_GTF                                                     } from '../../modules/local/filter_gtf/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_GTF                                   } from '../../modules/nf-core/icountmini/segment/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_FILTGTF                               } from '../../modules/nf-core/icountmini/segment/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED                     } from '../../modules/local/resolve_unannotated/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER         } from '../../modules/local/resolve_unannotated/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_REGIONS             } from '../../modules/local/resolve_unannotated/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS } from '../../modules/local/resolve_unannotated/main'

workflow PREPARE_GENOME {
    take:
    fasta                      // file: .fasta
    fasta_fai                  // file: .fai
    smrna_fasta                // file: .fasta
    smrna_fasta_fai            // file: .fai
    gtf                        // file: .gtf
    target_genome_index        // folder: index
    smrna_genome_index         // folder: index
    target_chrom_sizes         // file: .txt
    smrna_chrom_sizes          // file: .txt
    longest_transcript         // file: .txt
    longest_transcript_fai     // file: .fai
    longest_transcript_gtf     // file: .gtf
    filtered_gtf               // file: .gtf
    seg_gtf                    // file: .gtf
    seg_filt_gtf               // file: .gtf
    seg_resolved_gtf           // file: .gtf
    seg_resolved_gtf_genic     // file: .gtf
    regions_gtf                // file: .gtf
    regions_filt_gtf           // file: .gtf
    regions_resolved_gtf       // file: .gtf
    regions_resolved_gtf_genic // file: .gtf

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
    // MODULE: Uncompress genome smrna_fasta file if required
    //
    ch_smrna_fasta = Channel.empty()
    if (smrna_fasta.toString().endsWith(".gz")) {
        ch_smrna_fasta = GUNZIP_SMRNA_FASTA ( [ [id:smrna_fasta.baseName], smrna_fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_SMRNA_FASTA.out.versions)
    } else {
        ch_smrna_fasta = Channel.of([ [id:smrna_fasta.baseName], smrna_fasta ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fasta]
    //ch_smrna_fasta | view

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
    if (target_genome_index) {
        if (target_genome_index.toString().endsWith(".tar.gz")) {
            ch_star_index = UNTAR_STAR ( [ [:], target_genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_STAR.out.versions)
        } else {
            ch_star_index = Channel.of([ [:] , target_genome_index ])
        }
    }
    else {
        ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map{it[1]}, ch_gtf.map{it[1]} ).index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // MODULES: Uncompress Bowtie index or generate if required
    //
    ch_bt_index = Channel.empty()
    if (smrna_genome_index) {
        if (smrna_genome_index.toString().endsWith(".tar.gz")) {
            ch_bt_index = UNTAR_BT ( [ [:], smrna_genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BT.out.versions)
        } else {
            ch_bt_index = Channel.of([ [:] , smrna_genome_index ])
        }
    }
    else {
        ch_bt_index = BOWTIE_BUILD ( ch_smrna_fasta.map{it[1]} ).index
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    }

    //
    // MODULE: Create fasta fai if required
    //
    ch_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_fasta_fai = Channel.of([ [id:fasta_fai.baseName], fasta_fai ])
    } else {
        TARGET_INDEX (
            ch_fasta,
            [[],[]]
        )
        ch_fasta_fai = TARGET_INDEX.out.fai
        ch_versions = ch_versions.mix(TARGET_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Create fasta fai if required for smrna genome
    //
    ch_smrna_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_smrna_fasta_fai = Channel.of([ [id:smrna_fasta_fai.baseName], fasta_fai ])
    } else {
        SMRNA_INDEX (
            ch_smrna_fasta,
            [[],[]]
        )
        ch_smrna_fasta_fai = SMRNA_INDEX.out.fai
        ch_versions = ch_versions.mix(SMRNA_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Calc target chrom sizes
    //
    ch_target_chrom_sizes = target_chrom_sizes
    if(!target_chrom_sizes) {
        ch_target_chrom_sizes = TARGET_CHROM_SIZE ( ch_fasta ).sizes
        ch_versions  = ch_versions.mix(TARGET_CHROM_SIZE.out.versions)
    }

    //
    // MODULE: Calc smrna chrom sizes
    //
    ch_smrna_chrom_sizes = smrna_chrom_sizes
    if(!smrna_chrom_sizes) {
        ch_smrna_chrom_sizes = SMRNA_CHROM_SIZE ( ch_smrna_fasta ).sizes
        ch_versions  = ch_versions.mix(SMRNA_CHROM_SIZE.out.versions)
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
    // MODULE: Find the longest transcript from the primary genome
    //
    ch_longest_transcript     = Channel.empty()
    ch_longest_transcript_fai = Channel.empty()
    ch_longest_transcript_gtf = Channel.empty()
    if (longest_transcript && longest_transcript_fai && longest_transcript_gtf){
        ch_longest_transcript     =  longest_transcript
        ch_longest_transcript_fai =  longest_transcript_fai
        ch_longest_transcript_gtf =  longest_transcript_gtf
    } else {
        FIND_LONGEST_TRANSCRIPT (
            ch_gtf
        )
        ch_longest_transcript     = FIND_LONGEST_TRANSCRIPT.out.longest_transcript
        ch_longest_transcript_fai = FIND_LONGEST_TRANSCRIPT.out.longest_transcript_fai
        ch_longest_transcript_gtf = FIND_LONGEST_TRANSCRIPT.out.longest_transcript_gtf
        ch_versions               = ch_versions.mix(FIND_LONGEST_TRANSCRIPT.out.versions)
    }

    //
    // MODULE: Filter the GTF file
    //
    ch_filt_gtf = Channel.of( [ [id:filtered_gtf.baseName], filtered_gtf ] )
    if (!filtered_gtf){
        CLIPSEQ_FILTER_GTF (
            ch_gtf
        )
        ch_filt_gtf = CLIPSEQ_FILTER_GTF.out.gtf
        ch_versions = ch_versions.mix(CLIPSEQ_FILTER_GTF.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    // CLIPSEQ_FILTER_GTF.out.gtf | view

    //
    // MODULE: Segment GTF file using icount
    //
    ch_seg_gtf     = Channel.of( [ [id:seg_gtf.baseName], seg_gtf ] )
    ch_regions_gtf = Channel.of( [ [id:regions_gtf.baseName], regions_gtf ] )
   if (!seg_gtf || !regions_gtf){
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
    ch_seg_filt_gtf     = Channel.of( [ [id:seg_filt_gtf.baseName], seg_filt_gtf ] )
    ch_regions_filt_gtf = Channel.of( [ [id:regions_filt_gtf.baseName], regions_filt_gtf ] )
    if (!seg_filt_gtf || !regions_filt_gtf) {
        ICOUNT_SEG_FILTGTF (
            ch_filt_gtf,
            ch_fasta_fai.map{ it[1] }
        )
        ch_seg_filt_gtf     = ICOUNT_SEG_FILTGTF.out.gtf
        ch_regions_filt_gtf = ICOUNT_SEG_FILTGTF.out.regions
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ICOUNT_SEG_FILTGTF.out.gtf | view

    //
    // MODULE: Resolve the GTF regions that iCount did not annotate
    //
    ch_seg_resolved_gtf = Channel.of( [ [id:seg_resolved_gtf.baseName], seg_resolved_gtf ] )
    if (!params.seg_resolved_gtf) {
        RESOLVE_UNANNOTATED (
            ch_seg_gtf.map{ it[1] },
            ch_seg_filt_gtf.map{ it[1] },
            ch_gtf.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            false
        )
        ch_seg_resolved_gtf = RESOLVE_UNANNOTATED.out.gtf
        ch_versions         = ch_versions.mix(RESOLVE_UNANNOTATED.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED.out.gtf | view

    //
    // MODULE: Resolve the GTF regions that iCount did not annotate REGIONS FILE
    //
    ch_regions_resolved_gtf = Channel.of( [ [id:regions_resolved_gtf.baseName], regions_resolved_gtf ] )
    if (!regions_resolved_gtf) {
        RESOLVE_UNANNOTATED_REGIONS (
            ch_regions_gtf.map{ it[1] },
            ch_regions_filt_gtf.map{ it[1] },
            ch_gtf.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            false
        )
        ch_regions_resolved_gtf = RESOLVE_UNANNOTATED_REGIONS.out.gtf
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED_REGIONS.out.gtf | view

    //
    // MODULE: Resolve the GTF regions that iCount did not annotate with genic_other flag
    //
    ch_seg_resolved_gtf_genic = Channel.of( [ [id:seg_resolved_gtf_genic.baseName], seg_resolved_gtf_genic ] )
    if (!params.seg_resolved_gtf_genic) {
        RESOLVE_UNANNOTATED_GENIC_OTHER (
            ch_seg_gtf.map{ it[1] },
            ch_seg_filt_gtf.map{ it[1] },
            ch_gtf.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            true
        )
        ch_seg_resolved_gtf_genic = RESOLVE_UNANNOTATED_GENIC_OTHER.out.gtf
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED_GENIC_OTHER.out.gtf | view

    //
    // MODULE: Resolve the GTF regions that iCount did not annotate with genic_other flag REGIONS FILE
    //
    ch_regions_resolved_gtf_genic = Channel.of( [ [id:regions_resolved_gtf_genic.baseName], regions_resolved_gtf_genic ] )
    if (!params.regions_resolved_gtf_genic) {
        RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS (
            ch_regions_gtf.map{ it[1] },
            ch_regions_filt_gtf.map{ it[1] },
            ch_gtf.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            true
        )
        ch_regions_resolved_gtf_genic = RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS.out.gtf
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS.out.gtf | view


    emit:
    fasta                      = ch_fasta                      // channel: [ val(meta), [ fasta ] ]
    fasta_fai                  = ch_fasta_fai                  // channel: [ val(meta), [ fai ] ]
    smrna_fasta                = ch_smrna_fasta                // channel: [ val(meta), [ fasta ] ]
    smrna_fasta_fai            = ch_smrna_fasta_fai            // channel: [ val(meta), [ fai ] ]
    genome_index               = ch_star_index                 // channel: [ val(meta), [ star_index ] ]
    smrna_index                = ch_bt_index                   // channel: [ val(meta), [ bt2_index ] ]
    chrom_sizes                = ch_target_chrom_sizes         // channel: [ val(meta), [ txt ] ]
    smrna_chrom_sizes          = ch_smrna_chrom_sizes          // channel: [ val(meta), [ txt ] ]
    gtf                        = ch_gtf                        // channel: [ val(meta), [ gtf ] ]
    longest_transcript         = ch_longest_transcript         // channel: [ val(meta), [ txt ] ]
    longest_transcript_fai     = ch_longest_transcript_fai     // channel: [ val(meta), [ fai ] ]
    longest_transcript_gtf     = ch_longest_transcript_gtf     // channel: [ val(meta), [ fai ] ]
    filtered_gtf               = ch_filt_gtf                   // channel: [ val(meta), [ gtf ] ]
    seg_gtf                    = ch_seg_gtf                    // channel: [ val(meta), [ gtf ] ]
    seg_filt_gtf               = ch_seg_filt_gtf               // channel: [ val(meta), [ gtf ] ]
    seg_resolved_gtf           = ch_seg_resolved_gtf           // channel: [ val(meta), [ gtf ] ]
    seg_resolved_gtf_genic     = ch_seg_resolved_gtf_genic     // channel: [ val(meta), [ gtf ] ]
    regions_gtf                = ch_regions_gtf                // channel: [ val(meta), [ gtf ] ]
    regions_filt_gtf           = ch_regions_filt_gtf           // channel: [ val(meta), [ gtf ] ]
    regions_resolved_gtf       = ch_regions_resolved_gtf       // channel: [ val(meta), [ gtf ] ]
    regions_resolved_gtf_genic = ch_regions_resolved_gtf_genic // channel: [ val(meta), [ gtf ] ]
    versions                   = ch_versions                   // channel: [ versions.yml ]
}
