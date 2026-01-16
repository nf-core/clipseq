
include { sortAnnotationForBindingSiteFinder } from '../../../modules/local/BindingSiteFinder/sortAnnotationForBindingSiteFinder/main'
include { defineBindingSites } from '../../../modules/local/BindingSiteFinder/DefineBindingSites/main'
include { bsfQC } from '../../../modules/local/BindingSiteFinder/bsfQC/main'

workflow BINDINGSITEFINDER {

    take:

    gtf_ch // channel: [ val(meta), [ gtf ] ]
    bw_ch  // channel: [ val(meta), [ bigwig ] ]
    peak_ch // channel: [ peaks ]

    main:

    ch_versions = Channel.empty()

    // Destructure outputs from sortAnnotationForBindingSiteFinder(gtf_ch)
    tuple_ch = sortAnnotationForBindingSiteFinder(gtf_ch)
    gns_ch = tuple_ch.gns_rds
    regions_ch = tuple_ch.regions_rds
    ch_versions = ch_versions.mix(sortAnnotationForBindingSiteFinder.out.versions.first())

    // Pass annotation outputs, bigwig files and peak file to sortAnnotationForBindingSiteFinder(bw_ch)
    bs_ch = defineBindingSites(bw_ch, gns_ch, regions_ch, peak_ch)
    ch_versions = ch_versions.mix(defineBindingSites.out.versions.first())

    // Extract only the rds output from bs_ch
    rds_ch = bs_ch.binding_sites_rds

    // Pass the rds output to bsfQC
    qc_ch = bsfQC(rds_ch)
    ch_versions = ch_versions.mix(bsfQC.out.versions.first())

    emit:
    csv      = defineBindingSites.out.binding_sites_csv      // channel: [ val(meta), [ csv ] ]
    rds      = defineBindingSites.out.binding_sites_rds      // channel: [ val(meta), [ rds ] ]
    html     = bsfQC.out.bindingSiteFinderQC                 // channel: [ val(meta), [ html ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

