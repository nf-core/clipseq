process BindingSiteFinderAnalysis {
    tag "$meta.id"
    label 'process_low'

    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://melinak/bindingsitefinder:latest':
        'melinak/bindingsitefinder:latest' }"

    input:
        tuple val(meta), path(bw_files_folder)
        tuple val(meta), path(peaks)
        tuple val(meta), path(anno_gns)
        tuple val(meta), path(anno_regions)

    output:
        tuple val(meta), path("*binding_sites.rds"), emit: binding_sites_rds
        tuple val(meta), path("*binding_sites.csv"), emit: binding_sites_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the bindingSiteFinder sersion number will be printed to the console by the R script

    // optional parameters:
    // --peak_score_global_cuttoff $params.peak_score_global_cuttoff \\
    // --bsWidth $params.bsWidth \\
    // --peak_score_genewise_cuttoff $params.peak_score_genewise_cuttoff \\
    // --minWidth $params.minWidth \\
    // --minCrosslinks $params.minCrosslinks \\
    // --minCLSites $params.minCLSites \\
    // --maxBsWidth $params.maxBsWidth \\
    // // --reproducibility_cutoff $params.reproducibility_cutoff \\
    // // --reproducibility_nReps $params.reproducibility_nReps \\
    // --method_gene_overlaps $params.method_gene_overlaps \\
    // --rule_gene_overlaps $params.rule_gene_overlaps \\
    // --method_region_overlaps $params.method_region_overlaps \\
    // --rule_region_overlaps $params.rule_region_overlaps \\
    // --match_score $params.match_score \\
    // --match_geneID $params.match_geneID \\
    // --match_geneName $params.match_geneName \\
    // --match_geneType $params.match_geneType \\
    // --match_score_option $params.match_score_option \\

    """
    Rscript DefineBindingSites.R \\
        --bw_files_folder $bw_files_folder \\
        --peaks $peaks \\
        --anno_genes $anno_gns \\
        --anno_regions $anno_regions \\
        --sample_sheet $sample_sheet \\
        --output_path . \\
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch binding_sites.rds
    touch binding_sites.csv
    """
}