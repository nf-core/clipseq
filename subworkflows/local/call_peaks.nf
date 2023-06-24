//
// TODO: Summary
//

//
// MODULES
//
include { CLIPPY as CLIPPY_GENOME            } from "../../modules/nf-core/clippy/main"
include { ICOUNTMINI_SIGXLS                  } from "../../modules/nf-core/icountmini/sigxls/main"
include { ICOUNTMINI_PEAKS                   } from "../../modules/nf-core/icountmini/peaks/main"
include { GUNZIP as GUNZIP_ICOUNTMINI_SIGXLS } from "../../modules/nf-core/gunzip/main"
include { GUNZIP as GUNZIP_PEAKS_SIGXLS      } from "../../modules/nf-core/gunzip/main"
include { PARACLU as PARACLU_GENOME          } from "../../modules/nf-core/paraclu/main"
// include { PURECLIP } from "../modules/nf-core/pureclip/main.nf"

workflow CALL_PEAKS {
    take:
    callers   // channel: [ val ]
    xl_bed    // [ val(meta), [ bed ] ]
    gtf       // [ val(meta), [ gtf ] ]
    fasta_fai // [ val(meta), [ fasta, fai ] ]

    main:
    ch_versions = Channel.empty()

    ch_clippy_genome_peaks = Channel.empty()
    ch_clippy_genome_summits = Channel.empty()
    if('clippy' in callers) {
        CLIPPY_GENOME (
            xl_bed,
            gtf.collect{ it[1] },
            fasta_fai.collect{ it[1] }
        )
        
        ch_clippy_genome_peaks   = CLIPPY_GENOME.out.peaks
        ch_clippy_genome_summits = CLIPPY_GENOME.out.summits
        ch_versions              = ch_versions.mix(CLIPPY_GENOME.out.versions)
    }

    emit:
    clippy_genome_peaks   = ch_clippy_genome_peaks   // channel: [ val(meta), [ bedgraph ] ]
    clippy_genome_summits = ch_clippy_genome_summits // channel: [ val(meta), [ bedgraph ] ]
    versions              = ch_versions              // channel: [ versions.yml ]
}




        // if('icount' in callers) {

        //     ICOUNTMINI_SIGXLS (
        //         ch_target_crosslink_bed,
        //         ch_seg_resolved_gtf.collect{ it[1]}
                
        //     )

        //     ch_versions                      = ch_versions.mix(ICOUNTMINI_SIGXLS.out.versions)
        //     ch_icountmini_sigxls_gz          = ICOUNTMINI_SIGXLS.out.sigxls
        //     ch_icountmini_scores_gz          = ICOUNTMINI_SIGXLS.out.scores

        //     // CHANNEL: Create combined channel of input crosslinks and sigxls
        //     ch_peaks_input = ch_target_crosslink_bed
        //         .map{ [ it[0].id, it[0], it[1] ] }
        //         .join( ICOUNTMINI_SIGXLS.out.sigxls.map{ [ it[0].id, it[0], it[1] ] } )
        //         .map { [ it[1], it[2], it[4] ] }
        //     //EXAMPLE CHANNEL STRUCT: [ [id:test], BED(crosslinks), BED(sigxls) ]

        //     ICOUNTMINI_PEAKS (
        //         ch_peaks_input
        //     )

        //     ch_versions                      = ch_versions.mix(ICOUNTMINI_PEAKS.out.versions)
        //     ch_icountmini_peaks_gz           = ICOUNTMINI_PEAKS.out.peaks

        //     GUNZIP_ICOUNTMINI_SIGXLS (

        //         ch_icountmini_sigxls_gz

        //     )

        //     ch_versions                      = ch_versions.mix(GUNZIP_ICOUNTMINI_SIGXLS.out.versions)
        //     ch_icountmini_sigxls             = GUNZIP_ICOUNTMINI_SIGXLS.out.gunzip

        //     GUNZIP_ICOUNTMINI_PEAKS (

        //         ch_icountmini_peaks_gz

        //     )

        //     ch_versions                      = ch_versions.mix(GUNZIP_ICOUNTMINI_PEAKS.out.versions)
        //     ch_icountmini_peaks              = GUNZIP_ICOUNTMINI_PEAKS.out.gunzip

        // }

        // ch_paraclu_mincluster = Channel.value(params.paraclu_genome_params)

        // if('paraclu' in callers) {

        //     PARACLU_GENOME (
        //         ch_target_crosslink_bed,
        //         ch_paraclu_mincluster.collect()
        //     )

        //     ch_versions                      = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS.out.versions)
        //     ch_paraclu_genome_peaks          = PARACLU_GENOME.out.bed

        // }

        // // if('pureclip' in callers) {

        // // }
