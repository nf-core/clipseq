process CLIPSEQ_RESOLVE_UNANNOTATED {
    tag "$segmentation"
    label "process_single"

    conda "bioconda::pybedtools=0.9.0 conda-forge::plumbum=1.8.0"
    container "quay.io/goodwright/mulled-v2-9617f1b1a927f74fecc0c8b26ec773df8a8593b7:78688d5e3c856e3fbba8a63a7740b414dc4c0c5a-0"

    input:
    path segmentation
    path filt_segmentation
    path annotation
    path fai
    val genic_other

    output:
    path "*.gtf"         , emit: gtf
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name   = task.process
    output         = task.ext.output ?: "${filt_segmentation.simpleName}_genicOther${genic_other}.resolved.gtf"
    template 'resolve_unannotated.py'
}
