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

// Check non-mandatory input parameters to see if the files exist if they have been specified
def checkPathParamList = [
    params.multiqc_config,
    params.fasta_fai,
    params.ncrna_fasta_fai,
    params.genome_index,
    params.ncrna_genome_index,
    params.genome_chrom_sizes,
    params.ncrna_chrom_sizes,
    params.representative_transcript,
    params.representative_transcript_fai,
    params.representative_transcript_gtf,
    params.filtered_gtf,
    params.seg_gtf,
    params.regions_gtf,
    params.regions_filt_gtf,
    params.regions_resolved_gtf
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Define peak callers and check in list
caller_list = [ 'icount', 'paraclu', 'pureclip', 'clippy']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : []
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant caller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

// Creating a channel of callers to pass to CLIPQC later
ch_callers = Channel.value(callers)

// // Stage dummy files to be used as an optional input where required
ch_dummy_file  = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
ch_dummy_file2 = file("$projectDir/assets/dummy_file2.txt", checkIfExists: true)

// deseq2 header files:
ch_multiqc_merged_replicate_deseq2_pca_header        = file("$projectDir/assets/merged_replicate_deseq2_pca_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_deseq2_clustering_header = file("$projectDir/assets/merged_replicate_deseq2_clustering_header.txt", checkIfExists: true)

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

include { DUMP_SOFTWARE_VERSIONS                                          } from '../modules/local/dump_software_versions'
include { CLIPQC                                                          } from '../modules/local/clipqc'
include { LINUX_COMMAND as CONSENSUS_CROSSLINKS_REORDER_BED               } from '../modules/local/linux_command'
include { MERGE_SUMMARY                                                   } from '../modules/local/merge_summary'
include { ENCODE_MOVEUMI                                                  } from '../modules/local/encode_moveumi/main'
include { PEKA_MQC_PLOTS                                                  } from '../modules/local/peka_mqc_plots'
include { ICOUNTMINI_SUMMARY_MQC                                          } from '../modules/local/icountmini_summary_mqc'
include { MULTIMAP_CLASS_BINNING as BIN_ncRNA                             } from '../modules/local/multimap_class_binning'
include { MULTIMAP_CLASS_BINNING as BIN_REGIONS                           } from '../modules/local/multimap_class_binning'
include { GET_INIT_ALIGNED_XLINKS                                         } from '../modules/local/get_init_aligned_xlinks'
include { GET_INIT_ALIGNED_XLINKS as GET_INIT_ALIGNED_XLINKS_ncRNA        } from '../modules/local/get_init_aligned_xlinks'
include { IDENTIFY_UNBINNED                                               } from '../modules/local/identify_unbinned'
include { COMBINE_BINS                                                    } from '../modules/local/combine_bins'
include { rDNA_SUMMARY                                                    } from '../modules/local/rDNA_summary'
include { tRNA_SUMMARY                                                    } from '../modules/local/tRNA_summary'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME                                                        } from '../subworkflows/local/prepare_genome'
include { INPUT_CHECK                                                           } from '../subworkflows/local/input_check'
include { RNA_ALIGN                                                             } from '../subworkflows/local/rna_align'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as GENOME_UNIQUE_DEDUP                 } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as GENOME_MULTI_DEDUP                  } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as NCRNA_DEDUP                         } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { BAM_DEDUP_SAMTOOLS_UMICOLLAPSE as NCRNA_K1_DEDUP                      } from '../subworkflows/local/bam_dedup_samtools_umicollapse'
include { RESOLVE_GROUPS_AND_CROSSLINKS as GENOME_RESOLVE_GROUPS_AND_CROSSLINKS } from '../subworkflows/local/resolve_groupings_and_crosslinks'
include { RESOLVE_GROUPS_AND_CROSSLINKS as NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS  } from '../subworkflows/local/resolve_groupings_and_crosslinks'
include { RESOLVE_GROUPS_AND_CROSSLINKS as TE_RESOLVE_GROUPS_AND_CROSSLINKS     } from '../subworkflows/local/resolve_groupings_and_crosslinks'
include { TRANSCRIPTOME_PROCESSING                                              } from '../subworkflows/local/transcriptome_processing'
include { CONSENSUS_PEAK_TABLE as CLIPPY_CONSENSUS_PEAK_TABLE                   } from '../subworkflows/local/consensus_peak_table'
include { CONSENSUS_PEAK_TABLE as PARACLU_CONSENSUS_PEAK_TABLE                  } from '../subworkflows/local/consensus_peak_table'
include { CONSENSUS_PEAK_TABLE as ICOUNT_CONSENSUS_PEAK_TABLE                   } from '../subworkflows/local/consensus_peak_table'
include { MERGE_AND_SORT_TELESCOPE_BAMS                                         } from '../subworkflows/local/merge_and_sort_telescope_bams'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { BEDTOOLS_GROUPBY as CONSENSUS_CROSSLINKS_BEDTOOLS_GROUPBY } from '../modules/nf-core/bedtools/groupby/main'
include { LINUX_COMMAND as CONSENSUS_CROSSLINKS_SORT                } from '../modules/local/linux_command'
include { BEDTOOLS_MAP as CLIPPY_CONSENSUS_MAP                      } from '../modules/nf-core/bedtools/map/main'
include { CAT_CAT as CONSENSUS_CROSSLINKS_CAT_CAT                   } from '../modules/nf-core/cat/cat/main'
include { MULTIQC                                                   } from '../modules/nf-core/multiqc/main'
include { CLIPPY as CLIPPY_GENOME                                   } from "../modules/nf-core/clippy/main"
include { CLIPPY as CLIPPY_GENOME_CONSENSUS                         } from "../modules/nf-core/clippy/main"

include { ICOUNTMINI_SIGXLS                                         } from "../modules/nf-core/icountmini/sigxls/main"
include { ICOUNTMINI_PEAKS                                          } from "../modules/nf-core/icountmini/peaks/main"
include { GUNZIP as GUNZIP_ICOUNTMINI_SIGXLS                        } from "../modules/nf-core/gunzip/main"
include { GUNZIP as GUNZIP_ICOUNTMINI_PEAKS                         } from "../modules/nf-core/gunzip/main"

include { PARACLU as PARACLU_GENOME                                 } from "../modules/nf-core/paraclu/main"
include { PARACLU as PARACLU_GENOME_CONSENSUS                       } from "../modules/nf-core/paraclu/main"

include { PURECLIP as PURECLIP_WITH_CONTROL                         } from '../modules/nf-core/pureclip/main'
include { PURECLIP as PURECLIP_NO_CONTROL                           } from '../modules/nf-core/pureclip/main'
include { PEKA as PEKA_ICOUNT                                       } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_CLIPPY                                       } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_PARACLU                                      } from '../modules/nf-core/peka/main'
include { PEKA as PEKA_PURECLIP                                     } from '../modules/nf-core/peka/main'
include { ICOUNTMINI_SUMMARY                                        } from '../modules/nf-core/icountmini/summary/main'
include { ICOUNTMINI_METAGENE                                       } from '../modules/nf-core/icountmini/metagene/main'

include { BEDTOOLS_INTERSECT                                        } from '../modules/nf-core/bedtools/intersect/main'
include { TETRANSCRIPTS                                             } from '../modules/nf-core/tetranscripts/main'
include { TELESCOPE_ASSIGN                                          } from '../modules/nf-core/telescope/assign/main'
include { SAMTOOLS_VIEW as FILTER_UNIQUE_MAP_UPDATED                } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as FILTER_UNIQUE_MAP_OTHER                  } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SORT_BAMS_FOR_TELE                       } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SORT_INIT_ALIGNED_XLINKS                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as INDEX_INIT_ALIGNED_XLINKS               } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SORT_INIT_ALIGNED_XLINKS_ncRNA           } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as INDEX_INIT_ALIGNED_XLINKS_ncRNA         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as COUNT_PREMAPPED_READS                    } from '../modules/nf-core/samtools/view/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE                          } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'

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

    // Prepare mandatory params
    ch_input       = file(params.input)
    ch_fasta       = file(params.fasta)
    ch_ncrna_fasta = file(params.ncrna_fasta)
    ch_gtf         = file(params.gtf)

    // Prepare non-mandatory params
    ch_fasta_fai                         = []
    ch_ncrna_fasta_fai                   = []
    ch_genome_index                      = []
    ch_ncrna_genome_index                = []
    ch_genome_chrom_sizes                = []
    ch_ncrna_chrom_sizes                 = []
    ch_representative_transcript         = []
    ch_representative_transcript_fai     = []
    ch_representative_transcript_gtf     = []
    ch_filtered_gtf                      = []
    ch_seg_gtf                           = []
    ch_regions_gtf                       = []
    ch_regions_filt_gtf                  = []
    ch_regions_resolved_gtf              = []
    if(params.fasta_fai) { ch_fasta_fai = file(params.fasta_fai) }
    if(params.ncrna_fasta_fai) { ch_ncrna_fasta_fai = file(params.ncrna_fasta_fai) }
    if(params.genome_index) { ch_genome_index = file(params.genome_index) }
    if(params.ncrna_genome_index) { ch_ncrna_genome_index = file(params.ncrna_genome_index) }
    if(params.genome_chrom_sizes) { ch_genome_chrom_sizes = file(params.genome_chrom_sizes) }
    if(params.ncrna_chrom_sizes) { ch_ncrna_chrom_sizes = file(params.ncrna_chrom_sizes) }
    if(params.representative_transcript) { ch_representative_transcript = file(params.representative_transcript) }
    if(params.representative_transcript_fai) { ch_representative_transcript_fai = file(params.representative_transcript_fai) }
    if(params.representative_transcript_gtf) { ch_representative_transcript_gtf = file(params.representative_transcript_gtf) }
    if(params.filtered_gtf) { ch_filtered_gtf = file(params.filtered_gtf) }
    if(params.seg_gtf) { ch_seg_gtf = file(params.seg_gtf) }
    if(params.regions_gtf) { ch_regions_gtf = file(params.regions_gtf) }
    if(params.regions_filt_gtf) { ch_regions_filt_gtf = file(params.regions_filt_gtf) }
    if(params.regions_resolved_gtf) { ch_regions_resolved_gtf = file(params.regions_resolved_gtf) }

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
            ch_representative_transcript,
            ch_representative_transcript_fai,
            ch_representative_transcript_gtf,
            ch_filtered_gtf,
            ch_seg_gtf,
            ch_regions_gtf,
            ch_regions_filt_gtf,
            ch_regions_resolved_gtf,
            params.skip_filter_gtf,
            params.skip_transcriptome
        )
        ch_versions                          = ch_versions.mix(PREPARE_GENOME.out.versions)
        ch_fasta                             = PREPARE_GENOME.out.fasta.collect()
        ch_fasta_fai                         = PREPARE_GENOME.out.fasta_fai.collect()
        ch_gtf                               = PREPARE_GENOME.out.gtf.collect()
        ch_filtered_gtf                      = PREPARE_GENOME.out.filtered_gtf.collect()
        ch_genome_chrom_sizes                = PREPARE_GENOME.out.chrom_sizes.collect()
        ch_ncrna_fasta                       = PREPARE_GENOME.out.ncrna_fasta.collect()
        ch_ncrna_fasta_fai                   = PREPARE_GENOME.out.ncrna_fasta_fai.collect()
        ch_ncrna_chrom_sizes                 = PREPARE_GENOME.out.ncrna_chrom_sizes.collect()
        ch_representative_transcript         = PREPARE_GENOME.out.representative_transcript.collect()
        ch_representative_transcript_fai     = PREPARE_GENOME.out.representative_transcript_fai.collect()
        ch_representative_transcript_gtf     = PREPARE_GENOME.out.representative_transcript_gtf.collect()
        ch_seg_gtf                           = PREPARE_GENOME.out.seg_gtf.collect()
        ch_regions_gtf                       = PREPARE_GENOME.out.regions_gtf.collect()
        ch_regions_filt_gtf                  = PREPARE_GENOME.out.regions_filt_gtf.collect()
        ch_regions_resolved_gtf              = PREPARE_GENOME.out.regions_resolved_gtf.collect()
        ch_genome_index                      = PREPARE_GENOME.out.genome_index.collect()
        ch_ncrna_genome_index                = PREPARE_GENOME.out.ncrna_index.collect()
    }


    //
    // SUBWORKFLOW: Read in samplesheet, validate, stage input files and merge replicates
    //
    ch_fastq = Channel.empty()
    if(params.run_input_check) {
        INPUT_CHECK (
            ch_input,
            params.source
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
        ch_fastq    = INPUT_CHECK.out.reads
    }

    //EXAMPLE CHANNEL STRUCT: [[sample_name:h3k27me3_R1, group_name:h3k27me3, input_name:input, single_end:true, fastq:dsgsgh.fq.gz], [FASTQ]]
    // ch_fastq | view
    if(params.encode_eclip){
        ENCODE_MOVEUMI (
            ch_fastq
        )
        ch_versions = ch_versions.mix(ENCODE_MOVEUMI.out.versions)
        ch_fastq    = ENCODE_MOVEUMI.out.reads
    }
    //
    // SUBWORKFLOW: Extract UMI, trim and run before and after fastqc
    //
    if(params.source == "fastq" & params.run_preprocessing) {
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
    if(params.source == "fastq" & params.run_alignment) {

        ch_fastq
            .view { "ALIGN FASTQ: $it" }
        ch_ncrna_genome_index
            .view { "ALIGN ncRNA: $it" }
        ch_genome_index
            .view { "ALIGN GENOME INDEX: $it" }
        ch_gtf
            .view { "ALIGN GTF: $it" }
        ch_fasta
            .view { "ALIGN FASTA: $it" }

        RNA_ALIGN (
            ch_fastq,
            ch_ncrna_genome_index,
            ch_genome_index,
            ch_gtf,
            ch_fasta,
            params.skip_transcriptome
        )
        ch_versions                = ch_versions.mix(RNA_ALIGN.out.versions)
        ch_ncrna_bam               = RNA_ALIGN.out.ncrna_bam
        ch_ncrna_bai               = RNA_ALIGN.out.ncrna_bai
        ch_ncrna_log               = RNA_ALIGN.out.ncrna_log
        ch_ncrna_k1_bam            = RNA_ALIGN.out.ncrna_k1_bam
        ch_ncrna_k1_bai            = RNA_ALIGN.out.ncrna_k1_bai
        ch_genome_log              = RNA_ALIGN.out.genome_log_final
        ch_genome_unique_bam       = RNA_ALIGN.out.genome_unique_bam
        ch_genome_unique_bai       = RNA_ALIGN.out.genome_unique_bai
        ch_genome_multi_bam        = RNA_ALIGN.out.genome_multi_bam
        ch_genome_multi_bai        = RNA_ALIGN.out.genome_multi_bai
        ch_transcript_unique_bam   = RNA_ALIGN.out.transcript_unique_bam
        ch_transcript_unique_bai   = RNA_ALIGN.out.transcript_unique_bai
    }

    ch_paraclu_mincluster = Channel.value(params.paraclu_minValue)
    if(!params.skip_transcriptome) {
        TRANSCRIPTOME_PROCESSING(
            ch_transcript_unique_bam,
            ch_transcript_unique_bai,
            ch_representative_transcript,
            ch_representative_transcript_gtf,
            ch_representative_transcript_fai,
            callers,
            ch_paraclu_mincluster
        )
        ch_versions                      = ch_versions.mix(TRANSCRIPTOME_PROCESSING.out.versions)
        ch_transcript_bam                = TRANSCRIPTOME_PROCESSING.out.transcript_dedupe_bam
        ch_transcript_bai                = TRANSCRIPTOME_PROCESSING.out.transcript_dedupe_bai
        ch_trans_crosslink_bed           = TRANSCRIPTOME_PROCESSING.out.crosslink_bed
        ch_clippy_transcriptome_peaks    = TRANSCRIPTOME_PROCESSING.out.clippy_peaks
        ch_paraclu_transcriptome_peaks   = TRANSCRIPTOME_PROCESSING.out.paraclu_peaks
    }

    // Filling "final use" channels for when deduplication isn't run

    ch_genome_multi_final_bam = ch_genome_multi_bam
    ch_genome_multi_final_bai = ch_genome_multi_bai
    ch_genome_multi_final_bam_bai = ch_genome_multi_final_bam.join(ch_genome_multi_final_bai, by: 0)

    ch_genome_unique_final_bam = ch_genome_unique_bam
    ch_genome_unique_final_bai = ch_genome_unique_bai

    // Filling umi log as empty channel if dedup didn't run
    ch_umi_log = []

    // DEDUPLICATION //
    if(params.source == "fastq" & params.run_dedup) {
        // PREPARE CHANNELS
        ch_genome_unique_bam_bai = ch_genome_unique_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_unique_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_genome_multi_bam_bai = ch_genome_multi_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_multi_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_ncrna_bam_bai = ch_ncrna_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_ncrna_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_ncrna_k1_bam_bai = ch_ncrna_k1_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_ncrna_k1_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        /*
        * SUBWORKFLOW: Run umi deduplication on alignments
        */
        GENOME_UNIQUE_DEDUP (
            ch_genome_unique_bam_bai
        )
        ch_versions   = ch_versions.mix(GENOME_UNIQUE_DEDUP.out.versions)
        ch_genome_unique_final_bam = GENOME_UNIQUE_DEDUP.out.bam
        ch_genome_unique_final_bai = GENOME_UNIQUE_DEDUP.out.bai

        // NB multi is multi and unique
        GENOME_MULTI_DEDUP (
            ch_genome_multi_bam_bai
        )
        ch_versions   = ch_versions.mix(GENOME_MULTI_DEDUP.out.versions)
        ch_genome_multi_final_bam = GENOME_MULTI_DEDUP.out.bam
        ch_genome_multi_final_bai = GENOME_MULTI_DEDUP.out.bai
        ch_genome_multi_final_bam_bai = ch_genome_multi_final_bam.join(ch_genome_multi_final_bai, by: 0)
        ch_umi_log    = GENOME_MULTI_DEDUP.out.umi_log

        NCRNA_DEDUP (
            ch_ncrna_bam_bai
        )
        ch_versions   = ch_versions.mix(NCRNA_DEDUP.out.versions)

        NCRNA_K1_DEDUP (
            ch_ncrna_k1_bam_bai
        )
        ch_versions     = ch_versions.mix(NCRNA_K1_DEDUP.out.versions)
        ch_ncrna_k1_bam = NCRNA_K1_DEDUP.out.bam
        ch_ncrna_k1_bai = NCRNA_K1_DEDUP.out.bai
        ch_umi_ncrna_log      = NCRNA_K1_DEDUP.out.umi_log
    }

    // telescope gtf now required
    ch_rpt_annotation_gtf = Channel.value([[id: 'te_annotations'], file(params.ch_rpt_annotation_gtf, checkIfExists: true)])

    // Overall read assignment //
    // 1: get single xlink coordinates only

    GET_INIT_ALIGNED_XLINKS(
        ch_genome_multi_final_bam_bai
    )

    SORT_INIT_ALIGNED_XLINKS(
        GET_INIT_ALIGNED_XLINKS.out.bam,
        [[],[]]
    )

    ch_init_aligned_xlinks = SORT_INIT_ALIGNED_XLINKS.out.bam
        .map { meta, bam -> bam }
        .collect()
        .map { bams -> [[id: 'all_samples'], bams] }

    INDEX_INIT_ALIGNED_XLINKS(
        SORT_INIT_ALIGNED_XLINKS.out.bam
    )

    // 2: Binning to assign to TEs/ncRNA/tRNA

    BIN_ncRNA(
        ch_init_aligned_xlinks,
        ch_rpt_annotation_gtf,
    )

    // 3: Identify reads which were not assigned from 1st binning

    IDENTIFY_UNBINNED(
        ch_init_aligned_xlinks,
        BIN_ncRNA.out.discarded_reads
    )

    ch_init_unbinned_xlinks = IDENTIFY_UNBINNED.out.bam
        .map { meta, bam -> bam }
        .collect()
        .map { bams -> [[id: 'all_samples'], bams] }

    // 4: Run binning again on unassigned reads, at regional level (e.g. 5UTR | CDS | intron)

    BIN_REGIONS(
        ch_init_unbinned_xlinks,
        ch_regions_resolved_gtf
    )

    // 5: Count pre-mapped reads

    COUNT_PREMAPPED_READS(
        ch_ncrna_k1_bam.map { meta, bam -> [meta, bam, []] },
        [[],[]],
        []
    )

    ch_premapped_counts = COUNT_PREMAPPED_READS.out.bam
        .map { meta, bam -> bam}
        .collect()

    // 5: Combine assignments from both bins

    COMBINE_BINS(
        BIN_ncRNA.out.counted_reads,
        BIN_REGIONS.out.counted_reads,
        ch_premapped_counts
    )


    // tRNA summary insert //

    if(params.run_tRNA_summary) {

        ch_for_intersection = SORT_INIT_ALIGNED_XLINKS.out.bam
            .combine(ch_rpt_annotation_gtf.map { meta_gtf, gtf -> gtf })
            .map { meta, bed, gtf -> [meta, bed, gtf] }

        BEDTOOLS_INTERSECT(
            ch_for_intersection,
            [[],[]],
            "bed"
        )

        tRNA_SUMMARY(
           BEDTOOLS_INTERSECT.out.intersect
        )
    }

    // TE quantifcation insert //

    if(params.run_te) {
        // Check rmsk GTF has been provided
        if (!params.tetranscripts_gtf) {
            error "ERROR: --tetranscripts_gtf is required when --run_te is specified"
        }
        // Check rmsk GTF has been provided
        if (!params.ch_rpt_annotation_gtf) {
            error "ERROR: --ch_rpt_annotation_gtf is required when --run_te is specified"
        }

        ch_tetranscripts_gtf = Channel.value([[id: 'te_annotations'], file(params.tetranscripts_gtf, checkIfExists: true)])

        // Faux control .bam as not using DESeq2 aspect of TEtranscripts
        ch_bam_c = Channel.fromPath('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.bam')
             .map { file -> [[id: 'control'], file] }

        if (params.run_dedup){
        ch_all_t_bams = GENOME_MULTI_DEDUP.out.bam
            .map { meta, bam -> bam }
            .collect()
            .map { bams -> [[id: 'all_samples'], bams] }
        }

        if(!params.run_dedup){
        ch_all_t_bams = RNA_ALIGN.out.genome_multi_bam
            .map { meta, bam -> bam }
            .collect()
            .map { bams -> [[id: 'all_samples'], bams] }
        }

        TETRANSCRIPTS(
            ch_init_aligned_xlinks, // tx bam
            ch_bam_c, // control bam
            ch_gtf, // genome GTF
            ch_tetranscripts_gtf // te GTF
        )

        SORT_BAMS_FOR_TELE(
            SORT_INIT_ALIGNED_XLINKS.out.bam,
            [[],[]]
        )

        TELESCOPE_ASSIGN(
            SORT_BAMS_FOR_TELE.out.bam,
            ch_rpt_annotation_gtf
        )

        MERGE_AND_SORT_TELESCOPE_BAMS(
            TELESCOPE_ASSIGN.out.updated_bam,
            TELESCOPE_ASSIGN.out.other_bam
        )

    ch_genome_unique_final_bam = MERGE_AND_SORT_TELESCOPE_BAMS.out.bam
    ch_genome_unique_final_bai = MERGE_AND_SORT_TELESCOPE_BAMS.out.bai

    }

    // Crosslinks for ncRNA
    // genome-aligned files get crosslinks for bucketing, os ncRNA files get crosslinks here
    // to balance

    // join .bam and .bai into one channel

    ch_ncrna_k1_bam_bai_joined = ch_ncrna_k1_bam.join(ch_ncrna_k1_bai, by: 0)

    GET_INIT_ALIGNED_XLINKS_ncRNA(
        ch_ncrna_k1_bam_bai_joined
    )

    SORT_INIT_ALIGNED_XLINKS_ncRNA(
        GET_INIT_ALIGNED_XLINKS_ncRNA.out.bam,
        [[],[]]
    )

    INDEX_INIT_ALIGNED_XLINKS_ncRNA(
        SORT_INIT_ALIGNED_XLINKS_ncRNA.out.bam
    )

    ch_ncrna_k1_bam =     SORT_INIT_ALIGNED_XLINKS_ncRNA.out.bam
    ch_ncrna_k1_bai =     INDEX_INIT_ALIGNED_XLINKS_ncRNA.out.bai


    //
    // RESOLVE GROUPS AND GET CROSSLINKS: At this point, if groups have been specified, then we need to merge corresponding BAM files
    //
    // ch_genome_bam.view { item -> println("Pre-branch: $item") }
    if(params.run_crosslinking) {
        GENOME_RESOLVE_GROUPS_AND_CROSSLINKS(
            ch_genome_unique_final_bam ,
            ch_genome_unique_final_bai ,
            ch_fasta,
            ch_fasta_fai
        )
        ch_versions                                  = ch_versions.mix(GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.versions)
        ch_genome_crosslink_group_resolved_bed       = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_group_resolved
        ch_genome_crosslink_INDIVIDUAL_bed           = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_INDIVIDUAL
        ch_genome_crosslink_INDIVIDUAL_HASGROUP_bed  = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_INDIVIDUAL_HASGROUP
        ch_genome_crosslink_GROUP_HASGROUP_bed       = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_GROUP_HASGROUP
        ch_genome_peakcalling_bam                    = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.genome_peakcalling_bam
        ch_genome_peakcalling_bai                    = GENOME_RESOLVE_GROUPS_AND_CROSSLINKS.out.genome_peakcalling_bai

        NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS(
            ch_ncrna_k1_bam,
            ch_ncrna_k1_bai,
            ch_ncrna_fasta,
            ch_ncrna_fasta_fai
        )
        ch_versions                                    = ch_versions.mix(NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.versions)
        ch_ncrna_k1_crosslink_group_resolved_bed       = NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_group_resolved
        ch_ncrna_k1_crosslink_INDIVIDUAL_bed           = NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_INDIVIDUAL
        ch_ncrna_k1_crosslink_INDIVIDUAL_HASGROUP_bed  = NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_INDIVIDUAL_HASGROUP
        ch_ncrna_k1_crosslink_GROUP_HASGROUP_bed       = NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.crosslink_GROUP_HASGROUP


        // rDNA SUMMARY OUTPUT

        if(params.run_rDNA_summary) {
            rDNA_SUMMARY(
                NCRNA_RESOLVE_GROUPS_AND_CROSSLINKS.out.ncRNA_CPMnorm_bedgraph
                .map{ it[1] }
                .collect()
            )
        }

        // If filtering of GTF by transcripts is enabled, use the filtered GTF and its resolved regions, if not use those made by iCount-Mini
        ch_regions_used = params.skip_filter_gtf ? ch_regions_gtf : ch_regions_resolved_gtf
        ch_gtf_used = params.skip_filter_gtf ? ch_gtf : ch_filtered_gtf

        ICOUNTMINI_SUMMARY (
            ch_genome_crosslink_group_resolved_bed,
            ch_regions_used.map{ it[1] }
        )

        ch_merged_summaries = ICOUNTMINI_SUMMARY.out.summary_type
            .join( ICOUNTMINI_SUMMARY.out.summary_subtype, by: [0])
            .join( ICOUNTMINI_SUMMARY.out.summary_gene, by: [0])
            .join( ch_ncrna_k1_crosslink_group_resolved_bed, by: [0])
            .map { meta, type, subtype, gene, bed ->
                [meta, [type, subtype, gene, bed]]
            }


        MERGE_SUMMARY (
            ch_merged_summaries
        )

        ch_versions = ch_versions.mix(MERGE_SUMMARY.out.versions)
        ch_merged_summary_type = MERGE_SUMMARY.out.summary_type_adjusted
        ch_merged_summary_subtype = MERGE_SUMMARY.out.summary_subtype_adjusted
        ch_merged_summary_gene = MERGE_SUMMARY.out.summary_gene_adjusted

        // ICOUNTMINI_MERGE reformatting for mqc

        ICOUNTMINI_SUMMARY_MQC (
             MERGE_SUMMARY.out.summary_type_adjusted   .map { meta, f -> f }.collect(),
             MERGE_SUMMARY.out.summary_subtype_adjusted.map { meta, f -> f }.collect()
        )

        ICOUNTMINI_METAGENE (
            ch_genome_crosslink_group_resolved_bed,
            ch_regions_used.map{ it[1] }
        )

        if(params.consensus_peak){
            // Merge all xls into one file
            // want to use indivdual crosslinking files
            ch_all_crosslinks = ch_genome_crosslink_INDIVIDUAL_HASGROUP_bed.mix(ch_genome_crosslink_INDIVIDUAL_bed)

            ch_all_crosslinks
                .collect { it[1] }
                .map { crosslinks -> [[id:"allXL"], crosslinks]}
                .set { ch_consensus_crosslinks_bed }

            //ch_consensus_crosslinks_bed.view { item -> "consensus crosslinks: $item" }

            CONSENSUS_CROSSLINKS_CAT_CAT(
                ch_consensus_crosslinks_bed
            )

            CONSENSUS_CROSSLINKS_SORT(
                CONSENSUS_CROSSLINKS_CAT_CAT.out.file_out,
                [],
                false
            )
            // sum the counts to remove repeat entries
            CONSENSUS_CROSSLINKS_BEDTOOLS_GROUPBY(
                CONSENSUS_CROSSLINKS_SORT.out.file,
                5
            )

            CONSENSUS_CROSSLINKS_REORDER_BED(
                CONSENSUS_CROSSLINKS_BEDTOOLS_GROUPBY.out.bed,
                [],
                false
            )

            ch_consensus_crosslinks_final_bed = CONSENSUS_CROSSLINKS_REORDER_BED.out.file

        }
    }

    //
    // SUBWORKFLOW: Run peakcalling on genome
    //
    ch_clippy_genome_peaks          = Channel.empty()
    ch_icountmini_sigxls_gz         = Channel.empty()
    ch_icountmini_peaks_gz          = Channel.empty()
    ch_icountmini_sigxls            = Channel.empty()
    ch_paraclu_genome_peaks         = Channel.empty()

    //
    // PEKA CHANNELS
    //
    ch_peka_icountmini_peaks        = Channel.empty()
    ch_peka_clippy_peaks            = Channel.empty()
    ch_peka_paraclu_peaks           = Channel.empty()
    ch_peka_pureclip_peaks          = Channel.empty()


    // Create empty channel for DESEQ2_QC and PEKA_MQC_PLOTS for conditional use
    // initialise before the if(params.run_peakcalling) block
    ch_deseq2_qc_plots = Channel.empty()
    ch_peka_mqc_plots = Channel.empty()

    if(params.run_peakcalling) {

        if('clippy' in callers) {

            CLIPPY_GENOME (
                ch_genome_crosslink_group_resolved_bed,
                ch_gtf_used.map{ it[1] },
                ch_fasta_fai.map{ it[1] }
            )
            ch_versions             = ch_versions.mix(CLIPPY_GENOME.out.versions)
            ch_clippy_genome_peaks  = CLIPPY_GENOME.out.peaks

            if(params.consensus_peak){
                CLIPPY_GENOME_CONSENSUS (
                    ch_consensus_crosslinks_final_bed,
                    ch_gtf_used.map{ it[1] },
                    ch_fasta_fai.map{ it[1] }
                )

                CLIPPY_CONSENSUS_PEAK_TABLE (
                    ch_all_crosslinks,
                    CLIPPY_GENOME_CONSENSUS.out.peaks,
                    ch_fasta,
                    ch_fasta_fai,
                    ch_regions_used,
                    "Clippy_Consensus_AllCounts.tsv",
                    ch_multiqc_merged_replicate_deseq2_pca_header,
                    ch_multiqc_merged_replicate_deseq2_clustering_header,
                    params.skip_deseq2_qc
                )
                ch_versions = ch_versions.mix(CLIPPY_CONSENSUS_PEAK_TABLE.out.versions)
                ch_deseq2_qc_plots = ch_deseq2_qc_plots.mix(CLIPPY_CONSENSUS_PEAK_TABLE.out.deseq2_qc_plots)
            }

            if(params.run_peka) {
                PEKA_CLIPPY (
                    ch_clippy_genome_peaks,
                    ch_genome_crosslink_group_resolved_bed,
                    ch_fasta.map{ it[1] },
                    ch_fasta_fai.map{ it[1] },
                    ch_regions_used.map{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_CLIPPY.out.versions)
                ch_peka_mqc_plots = ch_peka_mqc_plots.mix(PEKA_CLIPPY.out.distribution)
            }

        }

        if('icount' in callers) {

            ICOUNTMINI_SIGXLS (
                ch_genome_crosslink_group_resolved_bed,
                ch_seg_gtf.map{ it[1]}

            )

            ch_versions                      = ch_versions.mix(ICOUNTMINI_SIGXLS.out.versions)
            ch_icountmini_sigxls_gz          = ICOUNTMINI_SIGXLS.out.sigxls

            // CHANNEL: Create combined channel of input crosslinks and sigxls
            ch_peaks_input = ch_genome_crosslink_group_resolved_bed
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

            if(params.consensus_peak){
                CONSENSUS_ICOUNTMINI_SIGXLS (
                    ch_consensus_crosslinks_final_bed,
                    ch_seg_gtf.map{ it[1]}
                )
                // CHANNEL: Create combined channel of input crosslinks and sigxls
                ch_consensus_peaks_input = ch_consensus_crosslinks_final_bed
                    .map{ [ it[0].id, it[0], it[1] ] }
                    .join( CONSENSUS_ICOUNTMINI_SIGXLS.out.sigxls.map{ [ it[0].id, it[0], it[1] ] } )
                    .map { [ it[1], it[2], it[4] ] }
                //EXAMPLE CHANNEL STRUCT: [ [id:test], BED(crosslinks), BED(sigxls) ]
                CONSENSUS_ICOUNTMINI_PEAKS (
                    ch_consensus_peaks_input
                )
                ch_consensus_peaks = CONSENSUS_ICOUNTMINI_PEAKS.out.peaks
                CONSENSUS_GUNZIP_ICOUNTMINI_SIGXLS (
                    CONSENSUS_ICOUNTMINI_SIGXLS.out.sigxls
                )
                CONSENSUS_GUNZIP_ICOUNTMINI_PEAKS (
                    CONSENSUS_ICOUNTMINI_PEAKS.out.peaks
                )
                ICOUNT_CONSENSUS_PEAK_TABLE (
                    ch_all_crosslinks,
                    CONSENSUS_GUNZIP_ICOUNTMINI_PEAKS.out.gunzip,
                    ch_fasta_fai,
                    ch_regions_used,
                    "iCount-Mini_Consensus_AllCounts.tsv"
                )
                ch_versions = ch_versions.mix(ICOUNT_CONSENSUS_PEAK_TABLE.out.versions)
            }

            if(params.run_peka) {
                PEKA_ICOUNT (
                    ch_icountmini_peaks,
                    ch_genome_crosslink_group_resolved_bed,
                    ch_fasta.map{ it[1] },
                    ch_fasta_fai.map{ it[1] },
                    ch_regions_used.map{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_ICOUNT.out.versions)
                ch_peka_mqc_plots = ch_peka_mqc_plots.mix(PEKA_ICOUNT.out.distribution)
            }

        }



        if('paraclu' in callers) {

            PARACLU_GENOME (
                ch_genome_crosslink_group_resolved_bed,
                ch_paraclu_mincluster
            )

            ch_versions                      = ch_versions.mix(PARACLU_GENOME.out.versions)
            ch_paraclu_genome_peaks          = PARACLU_GENOME.out.bed

            if(params.consensus_peak){
                PARACLU_GENOME_CONSENSUS (
                    ch_consensus_crosslinks_final_bed,
                    ch_paraclu_mincluster
                )

                PARACLU_CONSENSUS_PEAK_TABLE (
                    ch_all_crosslinks,
                    PARACLU_GENOME_CONSENSUS.out.bed,
                    ch_fasta,
                    ch_fasta_fai,
                    ch_regions_used,
                    "Paraclu_Consensus_AllCounts.tsv",
                    ch_multiqc_merged_replicate_deseq2_pca_header,
                    ch_multiqc_merged_replicate_deseq2_clustering_header,
                    params.skip_deseq2_qc
                )
                ch_versions = ch_versions.mix(PARACLU_CONSENSUS_PEAK_TABLE.out.versions)
                ch_deseq2_qc_plots = ch_deseq2_qc_plots.mix(PARACLU_CONSENSUS_PEAK_TABLE.out.deseq2_qc_plots)
            }

            if(params.run_peka) {
                PEKA_PARACLU (
                    ch_paraclu_genome_peaks,
                    ch_genome_crosslink_group_resolved_bed,
                    ch_fasta.map{ it[1] },
                    ch_fasta_fai.map{ it[1] },
                    ch_regions_used.map{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_PARACLU.out.versions)
                ch_peka_mqc_plots = ch_peka_mqc_plots.mix(PEKA_PARACLU.out.distribution)
            }
        }

        if('pureclip' in callers) {
            // Combine ch_genome_peakcalling_bam and ch_genome_peakcalling_bai correctly
            ch_genome_peakcalling = ch_genome_peakcalling_bam.join(ch_genome_peakcalling_bai, by: 0)

            // Print initial channel contents for debugging
            // ch_genome_peakcalling.view { item -> "Initial ch_genome_peakcalling item: $item" }

            ch_genome_peakcalling
                .branch { meta, bam, bai ->
                    control:     meta.control
                    no_control: !meta.control
                }
                .set { result }

            ch_genome_peakcalling
                .map{ meta, bam, bai -> [meta.id, bam, bai]
                }.set{ ch_genome_peakcalling_withid }

            result.control
                .map{ meta, bam, bai -> [meta.control, bam, bai, meta]
                }.set{ ch_genome_peakcalling_withControlid }

            ch_temp_pureclip_input = ch_genome_peakcalling_withControlid.join(ch_genome_peakcalling_withid, by: 0)
            // Structure is now [ControlID, IPBam, IPBai, IPMeta, Controlbam, Controlbai]

            // Check structure is what we expect
            //ch_temp_pureclip_input.view { item -> "Initial pureclip input merge channel: $item" }

            // Run PURECLIP for samples with control
            // reminder of input channel structure:
            //    tuple val(meta), path(ipbam), path(controlbam)
            //    tuple val(meta), path(ipbai), path(controlbai)
            //    tuple val(meta2), path(genome_fasta)
            //    val input_control
            ch_temp_pureclip_input
                .map{ ControlID, IPBam, IPBai, IPMeta, Controlbam, Controlbai -> [IPMeta, IPBam, Controlbam ]}
                .set{ ch_pureclip_bams_withcontrol }

            ch_temp_pureclip_input
                .map{ ControlID, IPBam, IPBai, IPMeta, Controlbam, Controlbai -> [IPMeta, IPBai, Controlbai ]}
                .set{ ch_pureclip_bais_withcontrol }

            PURECLIP_WITH_CONTROL(
                ch_pureclip_bams_withcontrol,
                ch_pureclip_bais_withcontrol,
                ch_fasta,
                true
            )

            // Run PURECLIP for samples without control
            result.no_control
                .map{ meta, bam, bai -> [meta, bam, [] ]}
                .set{ ch_pureclip_bams_nocontrol }

            result.no_control
                .map{ meta, bam, bai -> [meta, bai, [] ]}
                .set{ ch_pureclip_bais_nocontrol }

            PURECLIP_NO_CONTROL(
                ch_pureclip_bams_nocontrol,
                ch_pureclip_bais_nocontrol,
                ch_fasta,
                false
            )

            // Outputs from PURECLIP processes collected
            ch_versions                        = ch_versions.mix(PURECLIP_WITH_CONTROL.out.versions)
            ch_versions                        = ch_versions.mix(PURECLIP_NO_CONTROL.out.versions)
            ch_pureclip_genome_crosslinks      = PURECLIP_WITH_CONTROL.out.crosslinks.mix(PURECLIP_NO_CONTROL.out.crosslinks)
            ch_pureclip_genome_peaks           = PURECLIP_WITH_CONTROL.out.peaks.mix(PURECLIP_NO_CONTROL.out.peaks)

            if(params.run_peka) {
                // Need to make sure crosslinks and peaks are matched up correctly
                // After all the mixing
                ch_pureclip_genome_peaks.join(ch_genome_crosslink_group_resolved_bed, by: 0)
                    .set{ temp_matched_channel }

                temp_matched_channel
                    .map{ meta, peaks, crosslinks -> [meta, peaks] }
                    .set{ ch_pureclip_genome_peaks_matched }

                temp_matched_channel
                    .map{ meta, peaks, crosslinks -> [meta, crosslinks] }
                    .set{ ch_genome_crosslink_bed_matched }

                PEKA_PURECLIP(
                    ch_pureclip_genome_peaks_matched,
                    ch_genome_crosslink_bed_matched,
                    ch_fasta.map{ it[1] },
                    ch_fasta_fai.map{ it[1] },
                    ch_regions_used.map{ it[1] }
                )
                ch_versions = ch_versions.mix(PEKA_PURECLIP.out.versions)
                ch_peka_mqc_plots = ch_peka_mqc_plots.mix(PEKA_PURECLIP.out.distribution)
            }
        }
    }

    if(params.run_reporting) {

       // MODULE: Run clipqc

        CLIPQC (
            ch_callers,
            ch_ncrna_log.collect{ it[1] },
            ch_genome_log.collect{ it[1] },
            ch_umi_log.collect{ it[1] },
            ch_genome_crosslink_group_resolved_bed.collect{ it[1] },
            'icount' in callers ? ch_icountmini_peaks.collect{ it[1] } : [],
            'paraclu' in callers ? ch_paraclu_genome_peaks.collect {it[1]} : [],
            'clippy' in callers ? ch_clippy_genome_peaks.collect {it[1]} : [],
            'pureclip' in callers ? ch_pureclip_genome_peaks.collect {it[1] } : [],
            ch_merged_summary_type.collect{it[1]},
            ch_merged_summary_subtype.collect{it[1]},
            ch_merged_summary_gene.collect{it[1]},

        )

        // MODULE: Run peka multiqc plots

        PEKA_MQC_PLOTS (
            ch_peka_mqc_plots.collect{ it[1] }.flatten().collect()
        )

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

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_ncrna_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_genome_log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CLIPQC.out.tsv.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_deseq2_qc_plots.flatten().collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEKA_MQC_PLOTS.out.plots.flatten().collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(COMBINE_BINS.out.csv.collect{it[1]}.ifEmpty([]))

        if(params.run_crosslinking) {
            ch_multiqc_files = ch_multiqc_files.mix(ICOUNTMINI_SUMMARY_MQC.out.subtype.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ICOUNTMINI_SUMMARY_MQC.out.type.collect().ifEmpty([]))
        }

        ch_tmp = ch_ncrna_log.collect{it[1]}

        ch_multiqc_custom_config = ch_multiqc_custom_config.mix(COMBINE_BINS.out.config.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()

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

