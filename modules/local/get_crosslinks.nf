process GET_CROSSLINKS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.xl.bed")               , emit: bed
    tuple val(meta), path("*.xl.bedgraph")          , emit: bedgraph
    tuple val(meta), path("*.xl.CPMnorm.bedgraph")  , emit: norm_bedgraph
    path  "versions.yml"                            ,emit: versions


    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    bedtools bamtobed -i $bam >  shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{printf "%s\\t%d\\t%d\\t.\\t%.0f\\t+\\n", \$1, \$2, \$2+1, \$3}' > pos.bed
    bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{printf "%s\\t%d\\t%d\\t.\\t%.0f\\t-\\n", \$1, \$2, \$2+1, \$3}' > neg.bed
    cat pos.bed neg.bed | sort -k1,1 -k2,2n -k3,3 -k6,6 > ${prefix}.xl.bed
    cat ${prefix}.xl.bed | awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' > ${prefix}.xl.bedgraph
    TOTAL_VARIABLE=`cat ${prefix}.xl.bed | awk \'BEGIN {total=0} {total=total+\$5} END {print total}\'`
    cat ${prefix}.xl.bed | awk -v total=\$TOTAL_VARIABLE \'{printf "%s\\t%i\\t%i\\t%s\\t%f\\t%s\\n", \$1, \$2, \$3, \$4, 1000000*\$5/total, \$6}\' | awk \'{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}\' | sort -k1,1 -k2,2n > ${prefix}.xl.CPMnorm.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linux: NOVERSION
        bedtools: `bedtools --version | head -n 1`
    END_VERSIONS
    """
}
