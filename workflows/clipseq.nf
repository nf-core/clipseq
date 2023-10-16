/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowClipseq.initialise(params, log)

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    input: params.input,
    fasta: params.fasta,
    ncrna_fasta: params.ncrna_fasta,
    gtf: params.gtf
]
for (param in check_param_list) {
    if (!param.value) {
        exit 1, "Required parameter not specified: ${param.key}"
    }
    else {
        file(param.value, checkIfExists: true)
    }
}

// Check non-manditory input parameters to see if the files exist if they have been specified
def checkPathParamList = [
    params.multiqc_config,
    params.fasta_fai,
    params.ncrna_fasta_fai,
    params.genome_index,
    params.ncrna_genome_index,
    params.genome_chrom_sizes,
    params.ncrna_chrom_sizes,
    params.longest_transcript,
    params.longest_transcript_fai,
    params.longest_transcript_gtf,
    params.filtered_gtf,
    params.seg_gtf,
    params.seg_filt_gtf,
    params.seg_resolved_gtf,
    params.seg_resolved_gtf_genic,
    params.regions_gtf,
    params.regions_filt_gtf,
    params.regions_resolved_gtf,
    params.regions_resolved_gtf_genic
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Define peak callers and check in list
caller_list = [ 'icount', 'paraclu', 'pureclip', 'clippy']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : []
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant caller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

// // Stage dummy files to be used as an optional input where required
ch_dummy_file  = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
ch_dummy_file2 = file("$projectDir/assets/dummy_file2.txt", checkIfExists: true)

// // Check if an AWS iGenome has been provided to use the appropriate version of STAR
// def is_aws_igenome = false
// if (params.fasta && params.gtf) {
//     if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
//         is_aws_igenome = true
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules
//
include { SAMTOOLS_SIMPLE_VIEW as FILTER_TRANSCRIPTS } from '../modules/local/samtools_simple_view'
include { DUMP_SOFTWARE_VERSIONS                     } from '../modules/local/dump_software_versions'
include { CLIPQC                                     } from '../modules/local/clipqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME                                 } from '../subworkflows/local/prepare_genome'
include { PARSE_FASTQ_INPUT                              } from '../subworkflows/local/parse_fastq_input'
include { RNA_ALIGN                                      } from '../subworkflows/local/rna_align'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as GENOME_DEDUP } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as TRANS_DEDUP  } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { CALC_CROSSLINKS as CALC_GENOME_CROSSLINKS      } from '../subworkflows/local/calc_crosslinks'
include { CALC_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS  } from '../subworkflows/local/calc_crosslinks'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILT_TRANSCRIPT   } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILT_TRANSCRIPT } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                                          } from '../modules/nf-core/multiqc/main'
include { CLIPPY as CLIPPY_GENOME                          } from "../modules/nf-core/clippy/main"
include { CLIPPY as CLIPPY_TRANSCRIPTOME                   } from "../modules/nf-core/clippy/main"
include { ICOUNTMINI_SIGXLS                                } from "../modules/nf-core/icountmini/sigxls/main"
include { ICOUNTMINI_PEAKS                                 } from "../modules/nf-core/icountmini/peaks/main"
include { GUNZIP as GUNZIP_ICOUNTMINI_SIGXLS               } from "../modules/nf-core/gunzip/main"
include { GUNZIP as GUNZIP_PEAKS_SIGXLS                    } from "../modules/nf-core/gunzip/main"
include { GUNZIP as GUNZIP_ICOUNTMINI_PEAKS                } from "../modules/nf-core/gunzip/main"
include { PARACLU as PARACLU_GENOME                        } from "../modules/nf-core/paraclu/main"
include { PARACLU as PARACLU_TRANSCRIPTOME                 } from "../modules/nf-core/paraclu/main"
include { PURECLIP                                         } from '../modules/nf-core/pureclip/main'
include { PEKA as PEKA_ICOUNT                              } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_CLIPPY                              } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_PARACLU                             } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_PURECLIP                            } from '../modules/nf-core/peka/main'
include { ICOUNTMINI_SUMMARY                               } from '../modules/nf-core/icountmini/summary/main'
include { ICOUNTMINI_METAGENE                              } from '../modules/nf-core/icountmini/metagene/main'


//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CLIPSEQ {

    // Init
    ch_versions = Channel.empty()

    // Prepare manditory params
    ch_input       = file(params.input)
    ch_fasta       = file(params.fasta)
    ch_ncrna_fasta = file(params.ncrna_fasta)
    ch_gtf         = file(params.gtf)

    // Prepare non-manditory params
    ch_fasta_fai                  = []
    ch_ncrna_fasta_fai            = []
    ch_genome_index               = []
    ch_ncrna_genome_index         = []
    ch_genome_chrom_sizes         = []
    ch_ncrna_chrom_sizes          = []
    ch_longest_transcript         = []
    ch_longest_transcript_fai     = []
    ch_longest_transcript_gtf     = []
    ch_filtered_gtf               = []
    ch_seg_gtf                    = []
    ch_seg_filt_gtf               = []
    ch_seg_resolved_gtf           = []
    ch_seg_resolved_gtf_genic     = []
    ch_regions_gtf                = []
    ch_regions_filt_gtf           = []
    ch_regions_resolved_gtf       = []
    ch_regions_resolved_gtf_genic = []
    if(params.fasta_fai) { ch_fasta_fai = file(params.fasta_fai) }
    if(params.ncrna_fasta_fai) { ch_ncrna_fasta_fai = file(params.ncrna_fasta_fai) }
    if(params.genome_index) { ch_genome_index = file(params.genome_index) }
    if(params.ncrna_genome_index) { ch_ncrna_genome_index = file(params.ncrna_genome_index) }
    if(params.genome_chrom_sizes) { ch_genome_chrom_sizes = file(params.genome_chrom_sizes) }
    if(params.ncrna_chrom_sizes) { ch_ncrna_chrom_sizes = file(params.ncrna_chrom_sizes) }
    if(params.longest_transcript) { ch_longest_transcript = file(params.longest_transcript) }
    if(params.longest_transcript_fai) { ch_longest_transcript_fai = file(params.longest_transcript_fai) }
    if(params.longest_transcript_gtf) { ch_longest_transcript_gtf = file(params.longest_transcript_gtf) }
    if(params.filtered_gtf) { ch_filtered_gtf = file(params.filtered_gtf) }
    if(params.seg_gtf) { ch_seg_gtf = file(params.seg_gtf) }
    if(params.seg_filt_gtf) { ch_seg_filt_gtf = file(params.seg_filt_gtf) }
    if(params.seg_resolved_gtf) { ch_seg_resolved_gtf = file(params.seg_resolved_gtf) }
    if(params.seg_resolved_gtf_genic) { ch_seg_resolved_gtf_genic= file(params.seg_resolved_gtf_genic) }
    if(params.regions_gtf) { ch_regions_gtf = file(params.regions_gtf) }
    if(params.regions_filt_gtf) { ch_regions_filt_gtf = file(params.regions_filt_gtf) }
    if(params.regions_resolved_gtf) { ch_regions_resolved_gtf = file(params.regions_resolved_gtf) }
    if(params.regions_resolved_gtf_genic) { ch_regions_resolved_gtf_genic = file(params.regions_resolved_gtf_genic) }

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    if(params.run_genome_prep) {
        PREPARE_GENOME (
            ch_fasta,
            ch_fasta_fai,
            ch_ncrna_fasta,
            ch_ncrna_fasta_fai,
            ch_gtf,
            ch_genome_index,
            ch_ncrna_genome_index,
            ch_genome_chrom_sizes,
            ch_ncrna_chrom_sizes,
            ch_longest_transcript,
            ch_longest_transcript_fai,
            ch_longest_transcript_gtf,
            ch_filtered_gtf,
            ch_seg_gtf,
            ch_seg_filt_gtf,
            ch_seg_resolved_gtf,
            ch_seg_resolved_gtf_genic,
            ch_regions_gtf,
            ch_regions_filt_gtf,
            ch_regions_resolved_gtf,
            ch_regions_resolved_gtf_genic
        )
        ch_versions                   = ch_versions.mix(PREPARE_GENOME.out.versions)
        ch_fasta                      = PREPARE_GENOME.out.fasta
        ch_fasta_fai                  = PREPARE_GENOME.out.fasta_fai
        ch_gtf                        = PREPARE_GENOME.out.gtf
        ch_filtered_gtf               = PREPARE_GENOME.out.filtered_gtf
        ch_genome_chrom_sizes         = PREPARE_GENOME.out.chrom_sizes
        ch_ncrna_fasta                = PREPARE_GENOME.out.ncrna_fasta
        ch_ncrna_fasta_fai            = PREPARE_GENOME.out.ncrna_fasta_fai
        ch_ncrna_chrom_sizes          = PREPARE_GENOME.out.ncrna_chrom_sizes
        ch_longest_transcript         = PREPARE_GENOME.out.longest_transcript
        ch_longest_transcript_fai     = PREPARE_GENOME.out.longest_transcript_fai
        ch_longest_transcript_gtf     = PREPARE_GENOME.out.longest_transcript_gtf
        ch_seg_gtf                    = PREPARE_GENOME.out.seg_gtf
        ch_seg_filt_gtf               = PREPARE_GENOME.out.seg_filt_gtf
        ch_seg_resolved_gtf           = PREPARE_GENOME.out.seg_resolved_gtf
        ch_seg_resolved_gtf_genic     = PREPARE_GENOME.out.seg_resolved_gtf_genic
        ch_regions_gtf                = PREPARE_GENOME.out.regions_gtf
        ch_regions_filt_gtf           = PREPARE_GENOME.out.regions_filt_gtf
        ch_regions_resolved_gtf       = PREPARE_GENOME.out.regions_resolved_gtf
        ch_regions_resolved_gtf_genic = PREPARE_GENOME.out.regions_resolved_gtf_genic
        ch_genome_index        = PREPARE_GENOME.out.genome_index
        ch_ncrna_genome_index         = PREPARE_GENOME.out.ncrna_index
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate, stage input files and merge replicates
    //
    ch_fastq = Channel.empty()
    if(params.run_input_check) {
        PARSE_FASTQ_INPUT (
            ch_input
        )
        ch_versions = ch_versions.mix(PARSE_FASTQ_INPUT.out.versions)
        ch_fastq    = PARSE_FASTQ_INPUT.out.fastq
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:true], [FASTQ]]
    //ch_fastq | view

    //
    // SUBWORKFLOW: Extract UMI, trim and run b4 and after fastqc
    //
    if(params.run_preprocessing) {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_fastq,
            params.skip_fastqc,
            true,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
        ch_fastq    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:true], [FASTQ]]
    //ch_fastq | view

    //
    // SUBWORKFLOW: Align reads to ncrna and primary genomes
    //
    ch_ncrna_bam       = Channel.empty()
    ch_ncrna_bai       = Channel.empty()
    ch_ncrna_log       = Channel.empty()
    ch_genome_log      = Channel.empty()
    ch_genome_bam      = Channel.empty()
    ch_genome_bai      = Channel.empty()
    ch_genome_stats    = Channel.empty()
    ch_genome_flagstat = Channel.empty()
    ch_genome_idxstats = Channel.empty()
    ch_transcript_bam  = Channel.empty()
    ch_transcript_bai  = Channel.empty()
    if(params.run_alignment) {
        RNA_ALIGN (
            ch_fastq,
            ch_ncrna_genome_index,
            ch_genome_index,
            ch_filtered_gtf,
            ch_fasta
        )
        ch_versions         = ch_versions.mix(RNA_ALIGN.out.versions)
        ch_ncrna_bam        = RNA_ALIGN.out.ncrna_bam
        ch_ncrna_bai        = RNA_ALIGN.out.ncrna_bai
        ch_ncrna_log        = RNA_ALIGN.out.ncrna_log
        ch_genome_log       = RNA_ALIGN.out.genome_log_final
        ch_genome_bam       = RNA_ALIGN.out.genome_bam
        ch_genome_bai       = RNA_ALIGN.out.genome_bai
        ch_genome_stats     = RNA_ALIGN.out.genome_stats
        ch_genome_flagstat  = RNA_ALIGN.out.genome_flagstat
        ch_genome_idxstats  = RNA_ALIGN.out.genome_idxstats
        ch_transcript_bam   = RNA_ALIGN.out.transcript_bam
        ch_transcript_bai   = RNA_ALIGN.out.transcript_bai
    }
    //ch_genome_bam | view

    if(params.run_filtering) {
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
            ch_longest_transcript
        )
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS.out.versions)

        //
        // SUBWORKFLOW: sort, index filtered bam
        //
        SAMTOOLS_SORT_FILT_TRANSCRIPT ( FILTER_TRANSCRIPTS.out.bam )
        ch_versions       = ch_versions.mix(SAMTOOLS_SORT_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bam = SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam

        SAMTOOLS_INDEX_FILT_TRANSCRIPT ( SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam )
        ch_versions       = ch_versions.mix(SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bai = SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.bai
    }

    ch_trans_stats    = Channel.empty()
    ch_trans_flagstat = Channel.empty()
    ch_trans_idxstats = Channel.empty()
    ch_genome_umi_log = Channel.empty()
    ch_trans_umi_log  = Channel.empty()
    if(params.run_dedup) {
        //
        // CHANNEL: Combine bam and bai files on id
        //
        ch_genome_bam_bai = ch_genome_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        //
        // SUBWORKFLOW: Run umi deduplication on genome-level alignments
        //
        GENOME_DEDUP (
            ch_genome_bam_bai,
            ch_fasta
        )
        ch_versions        = ch_versions.mix(GENOME_DEDUP.out.versions)
        ch_genome_bam      = GENOME_DEDUP.out.bam
        ch_genome_bai      = GENOME_DEDUP.out.bai
        ch_genome_stats    = GENOME_DEDUP.out.stats
        ch_genome_flagstat = GENOME_DEDUP.out.flagstat
        ch_genome_idxstats = GENOME_DEDUP.out.idxstats
        ch_genome_umi_log  = GENOME_DEDUP.out.umi_log

        //
        // SUBWORKFLOW: Run umi deduplication on transcript-level alignments
        //
        TRANS_DEDUP (
            ch_transcript_bam_bai,
            ch_fasta
        )
        ch_versions        = ch_versions.mix(TRANS_DEDUP.out.versions)
        ch_transcript_bam  = TRANS_DEDUP.out.bam
        ch_transcript_bai  = TRANS_DEDUP.out.bai
        ch_genome_stats    = TRANS_DEDUP.out.stats
        ch_genome_flagstat = TRANS_DEDUP.out.flagstat
        ch_genome_idxstats = TRANS_DEDUP.out.idxstats
        ch_trans_umi_log   = TRANS_DEDUP.out.umi_log
    }

    ch_genome_crosslink_bed           = Channel.empty()
    ch_genome_crosslink_coverage      = Channel.empty()
    ch_genome_crosslink_coverage_norm = Channel.empty()
    ch_trans_crosslink_bed            = Channel.empty()
    ch_trans_crosslink_coverage       = Channel.empty()
    ch_trans_crosslink_coverage_norm  = Channel.empty()
    if(params.run_crosslinking) {
        //
        // SUBWORKFLOW: Run crosslink calculation for  genome
        //
        CALC_GENOME_CROSSLINKS (
            ch_genome_bam,
            ch_fasta_fai
        )
        ch_versions                       = ch_versions.mix(CALC_GENOME_CROSSLINKS.out.versions)
        ch_genome_crosslink_bed           = CALC_GENOME_CROSSLINKS.out.bed
        ch_genome_crosslink_coverage      = CALC_GENOME_CROSSLINKS.out.coverage
        ch_genome_crosslink_coverage_norm = CALC_GENOME_CROSSLINKS.out.coverage_norm

        ICOUNTMINI_SUMMARY (
            ch_genome_crosslink_bed,
            ch_regions_resolved_gtf.collect{ it[1] }
        )

        ICOUNTMINI_METAGENE (
            ch_genome_crosslink_bed,
            ch_regions_resolved_gtf.collect{ it[1] }
        )

        //
        // SUBWORKFLOW: Run crosslink calculation for transcripts
        //
        CALC_TRANSCRIPT_CROSSLINKS (
            ch_transcript_bam,
            ch_longest_transcript_fai.map{ [[id:it.baseName], it] }
        )
        ch_versions                      = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS.out.versions)
        ch_trans_crosslink_bed           = CALC_TRANSCRIPT_CROSSLINKS.out.bed
        ch_trans_crosslink_coverage      = CALC_TRANSCRIPT_CROSSLINKS.out.coverage
        ch_trans_crosslink_coverage_norm = CALC_TRANSCRIPT_CROSSLINKS.out.coverage_norm
    }

    //
    // SUBWORKFLOW: Run peakcalling on genome
    //
    ch_clippy_genome_peaks          = Channel.empty()
    ch_clippy_transcriptome_peaks   = Channel.empty()
    ch_icountmini_sigxls_gz         = Channel.empty()
    ch_icountmini_peaks_gz          = Channel.empty()
    ch_icountmini_sigxls            = Channel.empty()
    ch_paraclu_genome_peaks         = Channel.empty()
    ch_paraclu_transcriptome_peaks  = Channel.empty()

    if(params.run_peakcalling) {

        if('clippy' in callers) {

            CLIPPY_GENOME (
                ch_genome_crosslink_bed,
                ch_filtered_gtf.collect{ it[1] },
                ch_fasta_fai.collect{ it[1] }
            )
            
            ch_clippy_genome_peaks           = CLIPPY_GENOME.out.peaks

            CLIPPY_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_longest_transcript_gtf,
                ch_longest_transcript_fai
            )
            
            ch_clippy_transcriptome_peaks    = CLIPPY_TRANSCRIPTOME.out.peaks
            ch_versions                      = ch_versions.mix(CLIPPY_GENOME.out.versions)

            if(params.run_peka) {
                PEKA_CLIPPY (
                    ch_clippy_genome_peaks,
                    ch_genome_crosslink_bed,
                    ch_fasta.collect{ it[1] },
                    ch_fasta_fai.collect{ it[1] },
                    ch_regions_resolved_gtf.collect{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_CLIPPY.out.versions)
            }

        }

        if('icount' in callers) {

            ICOUNTMINI_SIGXLS (
                ch_genome_crosslink_bed,
                ch_seg_resolved_gtf.collect{ it[1]}
                
            )

            ch_versions                      = ch_versions.mix(ICOUNTMINI_SIGXLS.out.versions)
            ch_icountmini_sigxls_gz          = ICOUNTMINI_SIGXLS.out.sigxls

            // CHANNEL: Create combined channel of input crosslinks and sigxls
            ch_peaks_input = ch_genome_crosslink_bed
                .map{ [ it[0].id, it[0], it[1] ] }
                .join( ICOUNTMINI_SIGXLS.out.sigxls.map{ [ it[0].id, it[0], it[1] ] } )
                .map { [ it[1], it[2], it[4] ] }
            //EXAMPLE CHANNEL STRUCT: [ [id:test], BED(crosslinks), BED(sigxls) ]

            ICOUNTMINI_PEAKS (
                ch_peaks_input
            )

            ch_versions                      = ch_versions.mix(ICOUNTMINI_PEAKS.out.versions)
            ch_icountmini_peaks_gz           = ICOUNTMINI_PEAKS.out.peaks

            GUNZIP_ICOUNTMINI_SIGXLS (

                ch_icountmini_sigxls_gz

            )

            ch_versions                      = ch_versions.mix(GUNZIP_ICOUNTMINI_SIGXLS.out.versions)
            ch_icountmini_sigxls             = GUNZIP_ICOUNTMINI_SIGXLS.out.gunzip

            GUNZIP_ICOUNTMINI_PEAKS (

                ch_icountmini_peaks_gz

            )

            ch_versions                      = ch_versions.mix(GUNZIP_ICOUNTMINI_PEAKS.out.versions)
            ch_icountmini_peaks              = GUNZIP_ICOUNTMINI_PEAKS.out.gunzip

            if(params.run_peka) {
                PEKA_ICOUNT (
                    ch_icountmini_peaks,
                    ch_genome_crosslink_bed,
                    ch_fasta.collect{ it[1] },
                    ch_fasta_fai.collect{ it[1] },
                    ch_regions_resolved_gtf.collect{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_ICOUNT.out.versions)
            }

        }

        ch_paraclu_mincluster = Channel.value(params.paraclu_minValue)

        if('paraclu' in callers) {

            PARACLU_GENOME (
                ch_genome_crosslink_bed,
                ch_paraclu_mincluster
            )

            ch_versions                      = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS.out.versions)
            ch_paraclu_genome_peaks          = PARACLU_GENOME.out.bed

            PARACLU_TRANSCRIPTOME (
                ch_trans_crosslink_bed,
                ch_paraclu_mincluster
            )
            ch_paraclu_transcriptome_peaks          = PARACLU_TRANSCRIPTOME.out.bed

            if(params.run_peka) {
                PEKA_PARACLU(
                    ch_paraclu_genome_peaks,
                    ch_genome_crosslink_bed,
                    ch_fasta.collect{ it[1] },
                    ch_fasta_fai.collect{ it[1] },
                    ch_regions_resolved_gtf.collect{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_PARACLU.out.versions)
            }
        }

        if('pureclip' in callers) {
            ch_pureclip_input_bam = ch_genome_bam.combine(Channel.of(ch_dummy_file))
            ch_pureclip_input_bai = ch_genome_bai.combine(Channel.of(ch_dummy_file2))

            PURECLIP( 
                ch_pureclip_input_bam,
                ch_pureclip_input_bai,
                ch_fasta,
                false
            )
        }

        ch_versions                        = ch_versions.mix(PURECLIP.out.versions)
        ch_pureclip_genome_crosslinks      = PURECLIP.out.crosslinks
        ch_pureclip_genome_peaks           = PURECLIP.out.peaks

        if(params.run_peka) {
            PEKA_PURECLIP(
                ch_pureclip_genome_peaks,
                ch_genome_crosslink_bed,
                ch_fasta.collect{ it[1] },
                ch_fasta_fai.collect{ it[1] },
                ch_regions_resolved_gtf.collect{ it[1] }
            )
            ch_versions = ch_versions.mix(PEKA_PURECLIP.out.versions)
        }

    }

    if(params.run_reporting) {
        //
        // MODULE: Collect software versions
        //
        DUMP_SOFTWARE_VERSIONS (
            ch_versions.unique().collectFile()
        )

        //
        // MODULE: Run clipqc
        //
        // CLIPSEQ_CLIPQC (
        //     ch_bt_log.map{ it[1] },
        //     ch_star_log.map{ it[1] },
        //     ch_umi_log.map{ it[1] },
        //     ch_genome_crosslink_bed.map{ it[1] },
        //     ICOUNT_ANALYSE.out.bed_peaks.map{ it[1] },
        //     PARACLU_ANALYSE_GENOME.out.peaks.map{ it[1] },
        //     CLIPPY_GENOME.out.peaks.map{ it[1] }
        // )

        //
        // MODULE: Run multiqc
        //
        workflow_summary    = WorkflowClipseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowClipseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect())
        
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_ncrna_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_genome_log.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()

        // MULTIQC (
        //     ch_multiqc_config,
        //     DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect(),
        //     DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect(),
        //     ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yml"),
        //     FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
        //     FASTQC_TRIMGALORE.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]),
        //     FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
        //     ch_bt_log.collect{it[1]}.ifEmpty([]),
        //     ch_star_log.collect{it[1]}.ifEmpty([]),
        //     CLIPSEQ_CLIPQC.out.tsv.collect().ifEmpty([])
        // )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
