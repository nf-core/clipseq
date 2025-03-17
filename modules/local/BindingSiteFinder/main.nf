process BindingSiteFinderAnalysis {
    container 'melinak/bindingsitefinder:latest'
    
    input:
        path bw_files_folder
        path peaks
        path genome_annotation
        path sample_sheet

    output:
        path "binding_sites.rds"
        path "binding_sites.csv"

    script:
    """
    Rscript DefineBindingSites.R \\
        --bw_files_folder $bw_files_folder \\
        --peaks $peaks \\
        --anno_genes $anno_gns \\
        --anno_regions $anno_regions \\
        --sample_sheet $sample_sheet \\
        --output_path . \\
        --peak_score_global_cuttoff $params.peak_score_global_cuttoff \\
        --bsWidth $params.bsWidth \\
        --peak_score_genewise_cuttoff $params.peak_score_genewise_cuttoff \\
        --minWidth $params.minWidth \\
        --minCrosslinks $params.minCrosslinks \\
        --minCLSites $params.minCLSites \\
        --maxBsWidth $params.maxBsWidth \\
        --reproducibility_cutoff $params.reproducibility_cutoff \\
        --reproducibility_nReps $params.reproducibility_nReps \\
        --method_gene_overlaps $params.method_gene_overlaps \\
        --rule_gene_overlaps $params.rule_gene_overlaps \\
        --method_region_overlaps $params.method_region_overlaps \\
        --rule_region_overlaps $params.rule_region_overlaps \\
        --match_score $params.match_score \\
        --match_geneID $params.match_geneID \\
        --match_geneName $params.match_geneName \\
        --match_geneType $params.match_geneType \\
        --match_score_option $params.match_score_option \\
        binding_sites.rds \\
        binding_sites.csv
    """
}