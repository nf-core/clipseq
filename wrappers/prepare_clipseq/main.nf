#!/usr/bin/env nextflow

include { PREPARE_GENOME } from '../../subworkflows/local/prepare_genome.nf'

workflow  {
    // Prepare mandatory params
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
}
