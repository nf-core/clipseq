process CLIPSEQ_RESOLVE_UNANNOTATED {
    tag "$unfilt_regs"
    label "process_single"

    conda "bioconda::pybedtools=0.9.0 conda-forge::plumbum=1.8.0"
    container "quay.io/goodwright/mulled-v2-9617f1b1a927f74fecc0c8b26ec773df8a8593b7:78688d5e3c856e3fbba8a63a7740b414dc4c0c5a-0"

    input:
    tuple val(meta), path(unfilt_regs)
    tuple val(meta), path(filt_regs)
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path  "*.log"                 , emit: log
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name   = task.process
    output         = task.ext.output ?: "${filt_regs.simpleName}.resolved.gtf"
    template 'resolve_unannotated.py'
}
