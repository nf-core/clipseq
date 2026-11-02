process DESEQ2_QC {
    tag "$meta.id"
    label "process_single"

    container 'community.wave.seqera.io/library/bioconductor-deseq2_cairo_fontconfig_freetype_pruned:8f20d524c7811e46'

    input:
    tuple val(meta), path(counts)
    path ch_multiqc_merged_replicate_deseq2_pca_header
    path ch_multiqc_merged_replicate_deseq2_clustering_header

    output:
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*.rds"                , optional:true, emit: rds
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path "*.png"                , optional:true, emit: plots
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //template 'deseq2_qc.R'
    """
    Rscript $projectDir/modules/local/deseq2_qc/templates/deseq2_qc.R \
    -i ${counts} \
    -c 3 \
    -d "PeakID" \
    -o ${prefix} \
    --vst
    """

    stub:
    """
    touch ${meta.id}.pdf
    touch ${meta.id}.RData
    touch ${meta.id}.rds
    touch ${meta.id}.pca.vals.txt
    touch ${meta.id}.pca.vals_mqc.tsv
    touch ${meta.id}.sample.dists.txt
    touch ${meta.id}.sample.dists_mqc.tsv
    touch ${meta.id}.log
    path "size_factors"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """


}



