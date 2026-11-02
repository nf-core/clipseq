process PEKA_MQC_PLOTS {
    label "process_single"

    container 'community.wave.seqera.io/library/matplotlib_numpy_pandas_pathlib_pruned:2eddce7a3e9fc5bc'

    input:
    path ch_peka_mqc_plots

    output:
    path "*.png", optional: true, emit: plots
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''


    // Klara's PEKA MQC script
    """
    python ${projectDir}/modules/local/peka_mqc_plots/templates/peka_multiqc_summary.py "${ch_peka_mqc_plots.join(' ')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """

}
