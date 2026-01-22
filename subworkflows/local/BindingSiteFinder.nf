nextflow.enable.dsl = 2

include { sortAnnotationForBindingSiteFinder } from '../../modules/local/sortAnnotationForBindingSiteFinder/main.nf'
include { BindingSiteFinderAnalysis } from '../../modules/local/BindingSiteFinder/main.nf'

workflow BindingSiteFinder{
    take:
    gtf_ch
    // bw_files_folder
    // peaks


    main:
        sorted_ch = sortAnnotationForBindingSiteFinder(gtf_ch)
        
        // BindingSiteFinderAnalysis(
        //     bw_files_folder,
        //     peaks,
        //     sorted_ch.gns,
        //     sorted_ch.regions,
        //     sample_sheet,
        //     // Here start optional parameters
        //     peak_score_global_cuttoff,
        //     bsWidth,
        //     peak_score_genewise_cuttoff,
        //     minWidth,
        //     minCrosslinks,
        //     minCLSites,
        //     maxBsWidth,
        //     reproducibility_cutoff,
        //     reproducibility_nReps,
        //     method_gene_overlaps,
        //     rule_gene_overlaps,
        //     method_region_overlaps,
        //     rule_region_overlaps,
        //     match_score,
        //     match_geneID,
        //     match_geneName,
        //     match_geneType,
        //     match_score_option
        // )
}
