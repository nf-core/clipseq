process CLIPQC {
    label "process_single"

    conda (params.enable_conda ? "bioconda::pybedtools=0.9.0 pandas=1.5.0" : null)
    container "quay.io/goodwright/mulled-v2-9617f1b1a927f74fecc0c8b26ec773df8a8593b7:78688d5e3c856e3fbba8a63a7740b414dc4c0c5a-0"

    input:
    path("premap/*")
    path("mapped/*")
    path("collapse/*")
    path("xlinks/*")
    path("icount/*")
    path("paraclu/*")
    path("clippy/*")

    output:
    path "*.tsv"         , emit: tsv
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    template 'clipqc.py'
}
