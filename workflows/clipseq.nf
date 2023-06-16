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
    smrna_fasta: params.smrna_fasta,
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
    params.smrna_fasta_fai,
    params.target_genome_index,
    params.smrna_genome_index,
    params.target_chrom_sizes,
    params.smrna_chrom_sizes,
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


// // Check rRNA databases for sortmerna
// if (params.remove_ribo_rna) {
//     ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
//     if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
// }

// // Check if file with list of fastas is provided when running BBSplit
// if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
//     ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
//     if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
// }

// // Check alignment parameters
// def prepareToolIndices  = []
// if (!params.skip_bbsplit)   { prepareToolIndices << 'bbsplit'             }
// if (!params.skip_alignment) { prepareToolIndices << params.aligner        }
// if (params.pseudo_aligner)  { prepareToolIndices << params.pseudo_aligner }

// // Get RSeqC modules to run
// def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
// if (params.bam_csi_index) {
//     for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
//         if (rseqc_modules.contains(rseqc_module)) {
//             rseqc_modules.remove(rseqc_module)
//         }
//     }
// }

// // Stage dummy file to be used as an optional input where required
// ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

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

// ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
// ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { PARSE_FASTQ_INPUT } from '../subworkflows/local/parse_fastq_input'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { UMITOOLS_EXTRACT } from '../modules/nf-core/umitools/extract/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//


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
    ch_smrna_fasta = file(params.smrna_fasta)
    ch_gtf         = file(params.gtf)

    // Prepare non-manditory params
    ch_fasta_fai                  = []
    ch_smrna_fasta_fai            = []
    ch_target_genome_index        = []
    ch_smrna_genome_index         = []
    ch_target_chrom_sizes         = []
    ch_smrna_chrom_sizes          = []
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
    if(params.smrna_fasta_fai) { ch_smrna_fasta_fai = file(params.smrna_fasta_fai) }
    if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index) }
    if(params.smrna_genome_index) { ch_smrna_genome_index = file(params.smrna_genome_index) }
    if(params.target_chrom_sizes) { ch_target_chrom_sizes = file(params.target_chrom_sizes) }
    if(params.smrna_chrom_sizes) { ch_smrna_chrom_sizes = file(params.smrna_chrom_sizes) }
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
            ch_smrna_fasta,
            ch_smrna_fasta_fai,
            ch_gtf,
            ch_target_genome_index,
            ch_smrna_genome_index,
            ch_target_chrom_sizes,
            ch_smrna_chrom_sizes,
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
        ch_target_chrom_sizes         = PREPARE_GENOME.out.chrom_sizes
        ch_smrna_fasta                = PREPARE_GENOME.out.smrna_fasta
        ch_smrna_fasta_fai            = PREPARE_GENOME.out.smrna_fasta_fai
        ch_smrna_chrom_sizes          = PREPARE_GENOME.out.smrna_chrom_sizes
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
        ch_target_genome_index        = PREPARE_GENOME.out.genome_index
        ch_smrna_genome_index         = PREPARE_GENOME.out.smrna_index
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
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false], [FASTQ]]
    ch_fastq | view

    //
    // MODULE: Move umi to header if required
    //
    // if(params.run_move_umi_to_header){
    //     UMITOOLS_EXTRACT (
    //         ch_fastq
    //     )
    //     ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)
    //     ch_fastq    = UMITOOLS_EXTRACT.out.reads
    // }
    //UMITOOLS_EXTRACT.out.reads | view

    // if(params.run_trim_galore_fastqc) {
    //     /*
    //     * SUBWORKFLOW: Run fastqc and trimming
    //     */
    //     FASTQC_TRIMGALORE (
    //         ch_fastq,
    //         params.skip_fastqc,
    //         params.skip_trimming
    //     )
    //     ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    //     ch_fastq    = FASTQC_TRIMGALORE.out.fastq
    // }
    // //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false], [FASTQ]]
    // //ch_fastq | view


    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    // workflow_summary    = WorkflowClipseq.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowClipseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
