process ICOUNTMINI_SUMMARY_MQC{
     label 'process_single'

     container 'community.wave.seqera.io/library/r-tidyverse:2.0.0--386074aa810af1c3'

     input:
     path summary_types
     path summary_subtypes

     output:
     path "summary_subtype_qc.tsv", emit: subtype
     path "summary_type_qc.tsv", emit: type

     when:
     task.ext.when == null || task.ext.when

     script:
     template 'icount_mini_summary_reformat.R'

     stub:
     """
     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
     END_VERSIONS
     """
}
