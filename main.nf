#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/clipseq
========================================================================================
 nf-core/clipseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/clipseq
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/clipseq --input samplesheet.csv --fasta genome.fasta -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --      SET UP CONFIGURATION VARIABLE       -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Option for user supplied fasta and gtf and pipeline supplied smRNA
def smrna_list = ['human', 'mouse', 'rat', 'zebrafish', 'fruitfly', 'yeast']
if (!params.genome && params.smrna_org) {
    if (params.smrna_org in smrna_list) {
        params.smrna_fasta = params.smrna[ params.smrna_org ].smrna_fasta
    } else {
        params.smrna_fasta = false
        log.warn "There is no smRNA available for species '${params.smrna_org}'; pre-mapping will be skipped. Currently available options are: human, mouse, rat, fruitfly, zebrafish, yeast. Alternative you can supply your own smRNA fasta using --smrna_fasta"
    }
} else {
    if (params.genome && params.smrna.containsKey(params.genome)) {
        params.smrna_fasta = params.smrna[ params.genome ].smrna_fasta
    } else {
        params.smrna_fasta = false
        log.warn "There is no smRNA available for species '${params.genome}'; pre-mapping will be skipped. Currently available options are: human, mouse, rat, fruitfly, zebrafish, yeast. Alternative you can supply your own smRNA fasta using --smrna_fasta"
    }
}

// Auto-load genome files from genome config
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false

////////////////////////////////////////////////////
/* --    COLLECT CONFIGURATION PARAMETERS      -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input,
    params.fasta,
    params.gtf,
    params.star_index,
    params.smrna_fasta,
    params.fai
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if(!params.smrna_fasta) {
    if(params.genome) {
        log.warn "There is no available smRNA fasta file associated with the provided genome '${params.genome}'; pre-mapping will be skipped. A smRNA fasta file can be specified on the command line with --smrna_fasta"
    } else {
        log.warn "There is no smRNA fasta file suppled for genome specified; pre-mapping will be skipped. A smRNA fasta file can be specified on the command line with --smrna_fasta or --smrna_org"
    }
}

/*---- Check Peakcaller Options ---*/

callerList = [ 'icount', 'paraclu', 'pureclip', 'piranha']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : []
if ((callerList + callers).unique().size() != callerList.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${callerList.join(', ')}"
}

if ('icount' in callers) {
    icount_check = true
} else {
    icount_check = false
}

// Check genome is igenomes is used and icount peakcaller
// icount_compatible = [ 'GRCh37', 'GRCm38', 'TAIR10', 'EB2', 'UMD3.1', 'WBcel235', 'CanFam3.1', 'GRCz10', 'BDGP6', 'EquCab2', 'EB1', 'Galgal4', 'Gm01', 'Mmul_1', 'IRGSP-1.0', 'CHIMP2.1.4', 'Rnor_6.0', 'Rnor_5.0','R64-1-1', 'EF2', 'Sbi1', 'Sscrofa10.2', 'AGPv3' ]
icount_compatible = [] // Currently none of the iGenomes GTFs are compatible (even Ensembl - as different to the ones downloaded directly from Ensembl)
if (params.genome && ('icount' in callers) && !(params.genome in icount_compatible)) {
    icount_check = false
    log.warn "The provided genome '${params.genome}' is not compatible with the iCount peakcaller, so it will be skipped. Please see documentation"
}

// cannot run icount wihtout gtf file
if (!params.gtf && 'icount' in callers) {
    icount_check = false
    log.warn "iCount can only be run with a gtf annotation file - iCount will be skipped"
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

////////////////////////////////////////////////////
/* --               SET-UP INPUTS              -- */
////////////////////////////////////////////////////

params.umi_separator = params.move_umi ? '_' : ':' // Define default as ':' unless moving UMI with UMItools in which case '_'

if (params.smrna_fasta) ch_smrna_fasta = Channel.value(params.smrna_fasta)
if (params.star_index) ch_star_index = Channel.value(params.star_index)
if (params.gtf) ch_check_gtf = Channel.value(params.gtf)

// fai channels
if (params.fai) ch_fai_crosslinks = Channel.value(params.fai)
if (params.fai) ch_fai_icount = Channel.value(params.fai)
if (params.fai) ch_fai_icount_motif = Channel.value(params.fai)
if (params.fai) ch_fai_paraclu_motif = Channel.value(params.fai)
if (params.fai) ch_fai_size = Channel.value(params.fai)

// MultiQC empty channels from peakcaller checks
if (!('paraclu' in callers)) ch_paraclu_qc = Channel.empty()
if (!('icount' in callers) || !icount_check) ch_icount_qc = Channel.empty()
if (!('piranha' in callers)) ch_piranha_qc = Channel.empty()
if (!('pureclip' in callers)) ch_pureclip_qc = Channel.empty()

if (params.input) {
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true)
        .map{ row -> [ row.sample, file(row.fastq, checkIfExists: true) ] }
        .into{ ch_fastq; ch_fastq_fastqc_pretrim }
} else {
    exit 1, "Samples comma-separated input file not specified"
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --               HEADER LOG                 -- */
////////////////////////////////////////////////////

def summary = [:]
if (workflow.revision)                           summary['Pipeline Release'] = workflow.revision
summary['Run Name']                              = workflow.runName
summary['Input']                                 = params.input
if (params.fasta)                                summary['Fasta ref'] = params.fasta
if (params.gtf)                                  summary['GTF ref'] = params.gtf
if (params.star_index)                           summary['STAR index'] = params.star_index
if (params.save_index)                           summary['Save STAR index?'] = params.save_index
if (params.smrna_org)                            summary['SmallRNA organism ref'] = params.smrna_org
if (params.smrna_fasta)                          summary['SmallRNA ref'] = params.smrna_fasta
if (params.move_umi)                             summary['UMI pattern'] = params.move_umi
if (!params.skip_deduplication)                  summary['Deduplicate'] = params.deduplicate
if (!params.skip_deduplication && params.umi_separator) summary['UMI separator'] = params.umi_separator
if (params.peakcaller)                           summary['Peak caller'] = params.peakcaller
if (params.segment)                              summary['iCount segment'] = params.segment
if (icount_check)                                summary['Half window'] = params.half_window
if (icount_check)                                summary['Merge window'] = params.merge_window
if ('paraclu' in callers)                        summary['Min value'] = params.min_value
if ('paraclu' in callers)                        summary['Max density increase'] = params.min_density_increase
if ('paraclu' in callers)                        summary['Max cluster length'] = params.max_cluster_length
if ('pureclip' in callers)                       summary['Protein binding parameter'] = params.pureclip_bc
if ('pureclip' in callers)                       summary['Crosslink merge distance'] = params.pureclip_dm
if ('pureclip' in callers)                       summary['Chromosomes for HMM'] = params.pureclip_iv
if ('piranha' in callers)                        summary['Bin size'] = params.bin_size_both
if ('piranha' in callers)                        summary['Cluster distance'] = params.cluster_dist
summary['Max Resources']                         = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine)                    summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']                            = params.outdir
summary['Launch dir']                            = workflow.launchDir
summary['Working dir']                           = workflow.workDir
summary['Script dir']                            = workflow.projectDir
summary['User']                                  = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']                        = params.awsregion
    summary['AWS Queue']                         = params.awsqueue
    summary['AWS CLI']                           = params.awscli
}
summary['Config Profile']                        = workflow.profile
if (params.config_profile_description)           summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)               summary['Config Profile Contact'] = params.config_profile_contact
if (params.config_profile_url)                   summary['Config Profile URL'] = params.config_profile_url
summary['Config Files']                          = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']                    = params.email
    summary['E-mail on failure']                 = params.email_on_fail
    summary['MultiQC maxsize']                   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-clipseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/clipseq Workflow Summary'
    section_href: 'https://github.com/nf-core/clipseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    cutadapt --version > v_cutadapt.txt
    bowtie2 --version > v_bowtie2.txt
    STAR --version > v_star.txt
    samtools --version > v_samtools.txt
    umi_tools --version > v_umi_tools.txt
    bedtools --version > v_bedtools.txt
    preseq 2> v_preseq.txt
    # subread-align -v 2> v_subread.txt
    bam2fq.py --version > v_rseqc.txt
    iCount --version > v_icount.txt
    pureclip --version > v_pureclip.txt
    Piranha -about 2> v_piranha.txt
    echo "9" > v_paraclu.txt # Paraclu does not output a version
    meme -version > v_meme.txt
    python --version > v_python.txt
    pygmentize -V > v_pygments.txt
    pigz --version 2> v_pigz.txt
    perl -v > v_perl.txt

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

////////////////////////////////////////////////////
/* --             PREPROCESSING                -- */
////////////////////////////////////////////////////

/*
 * Generating premapping index
 */
if (params.smrna_fasta) {
    process generate_premap_index {
        tag "$smrna_fasta"
        label 'process_low'

        input:
        path(smrna_fasta) from ch_smrna_fasta

        output:
        path("${smrna_fasta.simpleName}.*.bt2") into ch_bt2_index

        script:
        """
        bowtie2-build --threads $task.cpus $smrna_fasta ${smrna_fasta.simpleName}
        """
    }
}

/*
 * Decompression
 */
// Recognise if fasta and/or gtf are compressed and decompress if so for STAR index generation
if (params.fasta) {
    if (hasExtension(params.fasta, 'gz')) {
        ch_fasta_gz = Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome reference fasta not found: ${params.fasta}" }
    } else {
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .into { ch_fasta; ch_fasta_fai; ch_fasta_dreme_icount; ch_fasta_dreme_paraclu; ch_fasta_pureclip; ch_fasta_dreme_pureclip; ch_fasta_dreme_piranha }
    }
}

if (params.fasta) {
    if (hasExtension(params.fasta, 'gz')) {
        process decompress_fasta {
            tag "$fasta_gz"
            label 'process_low'

            input:
            path(fasta_gz) from ch_fasta_gz

            output:
            path("*.fa") into (ch_fasta, ch_fasta_fai, ch_fasta_dreme_icount, ch_fasta_dreme_paraclu, ch_fasta_pureclip, ch_fasta_dreme_pureclip, ch_fasta_dreme_piranha)

            script:
            """
            pigz -d -c $fasta_gz > ${fasta_gz.baseName}
            """
        }
    }
}

/*
 * Generating fai index
 */
if (!params.fai) {
    process generate_fai {
        tag "$fasta"
        label 'process_low'

        input:
        path(fasta) from ch_fasta_fai

        output:
        path("*.fai") into (ch_fai_crosslinks, ch_fai_icount, ch_fai_icount_motif, ch_fai_paraclu_motif, ch_fai_pureclip_motif, ch_fai_piranha_motif)

        script:
        """
        samtools faidx $fasta
        """
    }
}

/*
 * Generating STAR index
 */
if (!params.star_index) {
    if (params.gtf) {
        if (hasExtension(params.gtf, 'gz')) {
            ch_gtf_gz_star = Channel
                .fromPath(params.gtf, checkIfExists: true)
                .ifEmpty { exit 1, "Genome reference gtf not found: ${params.gtf}" }
        } else {
            ch_gtf_star = Channel
                .fromPath(params.gtf, checkIfExists: true)
                .ifEmpty { exit 1, "Genome reference gtf not found: ${params.gtf}" }
        }
    }

    if (params.gtf) {
        if (hasExtension(params.gtf, 'gz')) {
            process decompress_gtf {
                tag "$gtf_gz"
                label 'process_low'

                input:
                path(gtf_gz) from ch_gtf_gz_star

                output:
                path("*.gtf") into ch_gtf_star

                script:
                """
                pigz -d -c $gtf_gz > ${gtf_gz.baseName}
                """
            }
        }
    }

    if (params.gtf) {
        process generate_star_index {
            tag "$fasta"
            label 'process_high'
            publishDir path: { params.save_index ? "${params.outdir}/STAR_index" : params.outdir },
                saveAs: { params.save_index ? it : null }, mode: params.publish_dir_mode

            input:
            path(fasta) from ch_fasta
            path(gtf) from ch_gtf_star

            output:
            path("STAR_${fasta.baseName}") into ch_star_index

            script:
            """
            samtools faidx $fasta
            NUM_BASES=`awk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

            mkdir STAR_${fasta.baseName}

            STAR \\
                --runMode genomeGenerate \\
                --runThreadN ${task.cpus} \\
                --genomeDir STAR_${fasta.baseName} \\
                --genomeFastaFiles $fasta \\
                --genomeSAindexNbases \$NUM_BASES \\
                --sjdbGTFfile $gtf
            """
        }
    } else if (!params.gtf) {
        process generate_star_index_no_gtf {
            tag "$fasta"
            label 'process_high'
            publishDir path: { params.save_index ? "${params.outdir}/STAR_index" : params.outdir },
                saveAs: { params.save_index ? it : null }, mode: params.publish_dir_mode

            input:
            path(fasta) from ch_fasta

            output:
            path("STAR_${fasta.baseName}") into ch_star_index

            script:
            """
            samtools faidx $fasta
            NUM_BASES=`awk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fasta}.fai`

            mkdir STAR_${fasta.baseName}

            STAR \\
                --runMode genomeGenerate --runThreadN ${task.cpus} \\
                --genomeDir STAR_${fasta.baseName} \\
                --genomeFastaFiles $fasta \\
                --genomeSAindexNbases \$NUM_BASES \\
            """
        }
    }
}


/*
 * Generating iCount segment file
 */
// iCount GTF input autodetects gz
if (params.peakcaller && icount_check) {
    if (!params.segment) {

        ch_gtf_icount = Channel
            .fromPath(params.gtf, checkIfExists: true)
            .ifEmpty { exit 1, "Genome reference gtf not found: ${params.gtf}" }

        process icount_segment {
            tag "$gtf"
            publishDir "${params.outdir}/icount", mode: params.publish_dir_mode

            input:
            path(gtf) from ch_gtf_icount
            path(fai) from ch_fai_icount

            output:
            path("icount_${gtf}") into ch_segment

            script:
            """
            mkdir tmp
            export ICOUNT_TMP_ROOT=\$PWD/tmp
            iCount segment $gtf icount_${gtf} $fai
            """
        }
    } else {
        ch_segment = Channel.value(params.segment)
    }
}

////////////////////////////////////////////////////
/* --             CLIP PIPELINE                -- */
////////////////////////////////////////////////////

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    tuple val(name), path(reads) from ch_fastq_fastqc_pretrim

    output:
    file "*fastqc.{zip,html}" into ch_fastqc_pretrim_mqc

    script:
    read_ext = reads.getName().split('\\.', 2)[1]
    read_name = reads.getName().split('\\.', 2)[0]
    new_reads = "${name}_reads_fastqc.${read_ext}"
    new_reads_simple = "${name}_reads_fastqc"
    """
    cp ${reads} ${new_reads}
    fastqc --quiet --threads $task.cpus ${new_reads}
    mv ${new_reads_simple}*.html ${name}_reads_fastqc.html
    mv ${new_reads_simple}*.zip ${name}_reads_fastqc.zip
    """
}

/*
 * STEP 1.1 - Move UMI to FastQ header if flagged
 */
if (params.move_umi) {
    process move_umi {
        tag "$name"
        label 'process_medium'
        publishDir "${params.outdir}/umi", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                    }

        input:
        tuple val(name), path(reads) from ch_fastq

        output:
        tuple val(name), path("${name}.umi.fastq.gz") into ch_umi_moved

        script:
        """
        umi_tools \\
            extract \\
            -p "$params.move_umi" \\
            -I $reads \\
            -S ${name}.umi.fastq.gz
        """
    }
} else {
    ch_umi_moved = ch_fastq
}

/*
 * STEP 2 - Read trimming
 */
process cutadapt {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/cutadapt", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads) from ch_umi_moved

    output:
    tuple val(name), path("${name}.trimmed.fastq.gz") into ch_trimmed
    path "*.log" into ch_cutadapt_mqc

    script:
    """
    [ -f ! ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
    cutadapt -j ${task.cpus} -a ${params.adapter} -m 12 -o ${name}.trimmed.fastq.gz ${name}.fastq.gz > ${name}_cutadapt.log
    """
}

/*
 * STEP 3 - Premapping
 */
if (params.smrna_fasta) {
    process premap {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/premap", mode: params.publish_dir_mode

        input:
        tuple val(name), path(reads) from ch_trimmed
        path(index) from ch_bt2_index.collect()

        output:
        tuple val(name), path("${name}.unmapped.fastq.gz") into ch_unmapped
        tuple val(name), path("${name}.premapped.bam"), path("${name}.premapped.bam.bai")
        path "*.log" into ch_premap_mqc, ch_premap_qc

        script:
        """
        bowtie2 -p $task.cpus -x ${index[0].simpleName} --un-gz ${name}.unmapped.fastq.gz -U $reads 2> ${name}.premap.log | \
        samtools sort -@ $task.cpus /dev/stdin > ${name}.premapped.bam && \
        samtools index -@ $task.cpus ${name}.premapped.bam
        """
    }
} else {
    ch_unmapped = ch_trimmed
    ch_premap_mqc = Channel.empty()
    ch_premap_qc = Channel.empty()
}

/*
 * STEP 4 - Aligning
 */
process align {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/mapped", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads) from ch_unmapped
    path(index) from ch_star_index.collect()

    output:
    tuple val(name), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai") into ch_aligned, ch_aligned_preseq
    path "*.Log.final.out" into ch_align_mqc, ch_align_qc

    script:
    clip_args = "--outFilterMultimapNmax 1 \
                --outFilterMultimapScoreRange 1 \
                --outSAMattributes All \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterType BySJout \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --outFilterScoreMin 10  \
                --alignEndsType Extend5pOfRead1 \
                --twopassMode Basic \
                --outSAMtype BAM Unsorted"
    """
    STAR \\
        --runThreadN $task.cpus \\
        --runMode alignReads \\
        --genomeDir $index \\
        --readFilesIn $reads --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${name}. $clip_args

    samtools sort -@ $task.cpus -o ${name}.Aligned.sortedByCoord.out.bam ${name}.Aligned.out.bam
    samtools index -@ $task.cpus ${name}.Aligned.sortedByCoord.out.bam
    """
}

/*
 * STEP 5 - Aligning QC
 */
process preseq {
    tag "$name"
    label 'process_low'
    publishDir "${params.outdir}/preseq", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam), path(bai) from ch_aligned_preseq

    output:
    path '*.ccurve.txt' into ch_preseq_mqc
    path '*.log'

    script:
    """
    preseq lc_extrap \\
        -output ${name}.ccurve.txt \\
        -verbose \\
        -bam \\
        -seed 42 \\
        $bam
    cp .command.err ${name}.command.log
    """
}

/*
 * STEP 6 - Deduplicate
 */
if (!params.skip_deduplication) {
    process dedup {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/dedup", mode: params.publish_dir_mode

        input:
        tuple val(name), path(bam), path(bai) from ch_aligned

        output:
        tuple val(name), path("${name}.dedup.bam"), path("${name}.dedup.bam.bai") into ch_dedup, ch_dedup_pureclip, ch_dedup_rseqc
        path "*.log" into ch_dedup_mqc, ch_dedup_qc

        script:
        """
        umi_tools \\
            dedup \\
            --umi-separator="$params.umi_separator" \\
            -I $bam \\
            -S ${name}.dedup.bam \\
            --output-stats=${name} \\
            --log=${name}.log
        samtools index -@ $task.cpus ${name}.dedup.bam
        """
    }
} else {
    ch_dedup = ch_aligned
    ch_dedup_mqc = Channel.empty()
    ch_dedup_qc = Channel.empty()
    ch_dedup_rseqc = ch_aligned
}

/*
 * STEP 6a - RSeQC
 */
if (params.gtf) {

    ch_gtf_rseqc = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .ifEmpty { exit 1, "Genome reference gtf not found: ${params.gtf}" }

    process rseqc {
        tag "$name"
        label 'process_low'
        publishDir "${params.outdir}/rseqc", mode: params.publish_dir_mode

        input:
        tuple val(name), path(bam), path(bai) from ch_dedup_rseqc
        path(gtf) from ch_gtf_rseqc.collect()

        output:
        path '*.read_distribution.txt' into ch_rseqc_mqc

        script:
        """
        gtf2bed $gtf > gene.bed

        read_distribution.py \\
            -i $bam \\
            -r gene.bed \\
            > ${name}.read_distribution.txt
        """
    }
} else {
    ch_rseqc_mqc = Channel.empty()
}

/*
 * STEP 7 - Identify crosslinks
 */
process get_crosslinks {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/xlinks", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam), path(bai) from ch_dedup
    path(fai) from ch_fai_crosslinks.collect()

    output:
    tuple val(name), path("${name}.xl.bed.gz") into ch_xlinks_icount, ch_xlinks_paraclu, ch_xlinks_piranha
    tuple val(name), path("${name}.xl.bedgraph.gz") into ch_xlinks_bedgraphs
    path "*.xl.bed.gz" into ch_xlinks_qc

    script:
    """
    bedtools bamtobed -i $bam > dedup.bed
    bedtools shift -m 1 -p -1 -i dedup.bed -g $fai > shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
    bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
    cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${name}.xl.bed.gz
    zcat ${name}.xl.bed.gz | awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' | pigz > ${name}.xl.bedgraph.gz
    """
}

/*
 * STEP 8a - Peak-call (iCount)
 */
if (params.peakcaller && icount_check) {
    process icount_peak_call {
        tag "$name"
        label 'process_low'
        publishDir "${params.outdir}/icount", mode: params.publish_dir_mode

        input:
        tuple val(name), path(xlinks) from ch_xlinks_icount
        path(segment) from ch_segment.collect()

        output:
        tuple val(name), path("${name}.${half_window}nt.sigxl.bed.gz") into ch_sigxls_icount
        tuple val(name), path("${name}.${half_window}nt_${merge_window}nt.peaks.bed.gz") into ch_peaks_icount
        path "*.peaks.bed.gz" into ch_icount_qc

        script:
        half_window = params.half_window
        merge_window = params.merge_window
        """
        mkdir tmp
        export ICOUNT_TMP_ROOT=\$PWD/tmp

        iCount peaks $segment $xlinks ${name}.${half_window}nt.sigxl.bed.gz --half_window ${half_window} --fdr 0.05

        pigz -d -c ${name}.${half_window}nt.sigxl.bed.gz | \\
        bedtools sort | \\
        bedtools merge -s -d ${merge_window} -c 4,5,6 -o distinct,sum,distinct | \\
        pigz > ${name}.${half_window}nt_${merge_window}nt.peaks.bed.gz
        """
    }

    if (params.motif) {
        process icount_motif_dreme {
            tag "$name"
            label 'process_low'
            publishDir "${params.outdir}/icount_motif", mode: params.publish_dir_mode

            input:
            tuple val(name), path(peaks) from ch_peaks_icount
            path(fasta) from ch_fasta_dreme_icount.collect()
            path(fai) from ch_fai_icount_motif.collect()

            output:
            tuple val(name), path("${name}_dreme/*") into ch_motif_dreme_icount

            script:
            motif_sample = params.motif_sample
            """
            pigz -d -c $peaks | awk '{OFS="\t"}{if(\$6 == "+") print \$1, \$2, \$2+1, \$4, \$5, \$6; else print \$1, \$3-1, \$3, \$4, \$5, \$6}' | \\
            bedtools slop -s -l 20 -r 20 -i /dev/stdin -g $fai | \\
            shuf -n $motif_sample > resized_peaks.bed

            bedtools getfasta -fi $fasta -bed resized_peaks.bed -fo resized_peaks.fasta

            dreme -norc -o ${name}_dreme -p resized_peaks.fasta
            """
        }
    }
}

/*
 * STEP 8b - Peak-call (paraclu)
 */
if ('paraclu' in callers) {
    process paraclu_peak_call {
        tag "$name"
        label 'process_low'
        publishDir "${params.outdir}/paraclu", mode: params.publish_dir_mode

        when:
        'paraclu' in callers

        input:
        tuple val(name), path(xlinks) from ch_xlinks_paraclu

        output:
        tuple val(name), path("${name}.${min_value}_${max_cluster_length}nt_${min_density_increase}.peaks.bed.gz") into ch_peaks_paraclu
        path "*.peaks.bed.gz" into ch_paraclu_qc

        script:
        min_value = params.min_value
        min_density_increase = params.min_density_increase
        max_cluster_length = params.max_cluster_length
        """
        pigz -d -c $xlinks | \\
        awk '{OFS = "\t"}{print \$1, \$6, \$3, \$5}' | \\
        sort -k1,1 -k2,2 -k3,3n > paraclu_input.tsv

        paraclu ${min_value} paraclu_input.tsv | \\
        paraclu-cut -d ${min_density_increase} -l ${max_cluster_length} | \\
        awk '{OFS = "\t"}{print \$1, \$3-1, \$4, ".", \$6, \$2}' | \\
        bedtools sort | \\
        pigz > ${name}.${min_value}_${max_cluster_length}nt_${min_density_increase}.peaks.bed.gz
        """
    }

    if (params.motif) {
        process paraclu_motif_dreme {
            tag "$name"
            label 'process_low'
            publishDir "${params.outdir}/paraclu_motif", mode: params.publish_dir_mode

            when:
            'paraclu' in callers

            input:
            tuple val(name), path(peaks) from ch_peaks_paraclu
            path(fasta) from ch_fasta_dreme_paraclu.collect()
            path(fai) from ch_fai_paraclu_motif.collect()

            output:
            tuple val(name), path("${name}_dreme/*") into ch_motif_dreme_paraclu

            script:
            motif_sample = params.motif_sample
            """
            pigz -d -c $peaks | awk '{OFS="\t"}{if(\$6 == "+") print \$1, \$2, \$2+1, \$4, \$5, \$6; else print \$1, \$3-1, \$3, \$4, \$5, \$6}' | \\
            bedtools slop -s -l 20 -r 20 -i /dev/stdin -g $fai | \\
            shuf -n $motif_sample > resized_peaks.bed

            bedtools getfasta -fi $fasta -bed resized_peaks.bed -fo resized_peaks.fasta

            dreme -norc -o ${name}_dreme -p resized_peaks.fasta
            """
        }
    }
}

/*
 * STEP 8c - Peak-call (PureCLIP)
 */
if ('pureclip' in callers) {
    process pureclip_peak_call {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/pureclip", mode: params.publish_dir_mode

        when:
        'pureclip' in callers

        input:
        tuple val(name), path(bam), path(bai) from ch_dedup_pureclip
        path(fasta) from ch_fasta_pureclip.collect()

        output:
        tuple val(name), path("${name}.sigxl.bed.gz") into ch_sigxlinks_pureclip
        tuple val(name), path("${name}.${dm}nt.peaks.bed.gz") into ch_peaks_pureclip
        path "*.peaks.bed.gz" into ch_pureclip_qc

        script:
        dm = params.pureclip_dm
        args = " -bc " + params.pureclip_bc
        args += " -dm " + params.pureclip_dm
        if (params.pureclip_iv) args += " -iv '" + params.pureclip_iv + "' "
        """
        pureclip \\
            -i $bam \\
            -bai $bai \\
            -g $fasta \\
            -nt $task.cpus \\
            $args \\
            -o "${name}.sigxl.bed" \\
            -or "${name}.${dm}nt.peaks.bed"

        pigz ${name}.sigxl.bed ${name}.${dm}nt.peaks.bed
        """
    }

    if (params.motif) {
        process pureclip_motif_dreme {
            tag "$name"
            label 'process_low'
            publishDir "${params.outdir}/pureclip_motif", mode: params.publish_dir_mode

            input:
            tuple val(name), path(peaks) from ch_peaks_pureclip
            path(fasta) from ch_fasta_dreme_pureclip.collect()
            path(fai) from ch_fai_pureclip_motif.collect()

            output:
            tuple val(name), path("${name}_dreme/*") into ch_motif_dreme_pureclip

            script:
            motif_sample = params.motif_sample
            """
            pigz -d -c $peaks | awk '{OFS="\t"}{if(\$6 == "+") print \$1, \$2, \$2+1, \$4, \$5, \$6; else print \$1, \$3-1, \$3, \$4, \$5, \$6}' | \\
            bedtools slop -s -l 20 -r 20 -i /dev/stdin -g $fai | \\
            shuf -n $motif_sample > resized_peaks.bed

            bedtools getfasta -fi $fasta -bed resized_peaks.bed -fo resized_peaks.fasta

            dreme -norc -o ${name}_dreme -p resized_peaks.fasta
            """
        }
    }
}


/*
 * STEP 8d - Peak-call (Piranha)
 */
if ('piranha' in callers) {
    process piranha_peak_call {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/piranha", mode: params.publish_dir_mode

        when:
        'piranha' in callers

        input:
        tuple val(name), path(xlinks) from ch_xlinks_piranha

        output:
        tuple val(name), path("${name}.${bin_size_both}nt_${cluster_dist}nt.peaks.bed.gz") into ch_peaks_piranha
        path "*.peaks.bed.gz" into ch_piranha_qc

        script:
        bin_size_both = params.bin_size_both
        cluster_dist = params.cluster_dist
        """
        pigz -d -c $xlinks | \\
        awk '{OFS="\t"}{for(i=0;i<\$5;i++) print }' \\
        > expanded.bed

        Piranha \\
            expanded.bed \\
            -s \\
            -b $bin_size_both \\
            -u $cluster_dist \\
            -o paraclu.bed

        awk '{OFS="\t"}{print \$1, \$2, \$3, ".", \$5, \$6}' paraclu.bed | \\
        pigz > ${name}.${bin_size_both}nt_${cluster_dist}nt.peaks.bed.gz
        """
    }

    if (params.motif) {
        process piranha_motif_dreme {
            tag "$name"
            label 'process_low'
            publishDir "${params.outdir}/piranha_motif", mode: params.publish_dir_mode

            input:
            tuple val(name), path(peaks) from ch_peaks_piranha
            path(fasta) from ch_fasta_dreme_piranha.collect()
            path(fai) from ch_fai_piranha_motif.collect()

            output:
            tuple val(name), path("${name}_dreme/*") into ch_motif_dreme_piranha

            script:
            motif_sample = params.motif_sample
            """
            pigz -d -c $peaks | awk '{OFS="\t"}{if(\$6 == "+") print \$1, \$2, \$2+1, \$4, \$5, \$6; else print \$1, \$3-1, \$3, \$4, \$5, \$6}' | \\
            bedtools slop -s -l 20 -r 20 -i /dev/stdin -g $fai | \\
            shuf -n $motif_sample > resized_peaks.bed

            bedtools getfasta -fi $fasta -bed resized_peaks.bed -fo resized_peaks.fasta

            dreme -norc -o ${name}_dreme -p resized_peaks.fasta
            """
        }
    }
}

/*
 * STEP 8 - QC plots
 */
process clipqc {
    label 'process_low'
    publishDir "${params.outdir}/clipqc", mode: params.publish_dir_mode

    input:
    file ('premap/*') from ch_premap_qc.collect().ifEmpty([])
    file ('mapped/*') from ch_align_qc.collect().ifEmpty([])
    file ('dedup/*') from ch_dedup_qc.collect().ifEmpty([])
    file ('xlinks/*') from ch_xlinks_qc.collect().ifEmpty([])
    path ('icount/*') from ch_icount_qc.collect().ifEmpty([])
    file ('paraclu/*') from ch_paraclu_qc.collect().ifEmpty([])
    file ('pureclip/*') from ch_pureclip_qc.collect().ifEmpty([])
    file ('piranha/*') from ch_piranha_qc.collect().ifEmpty([])

    output:
    path "*.tsv" into ch_clipqc_mqc

    script:
    clip_qc_args = ''

    if ('icount' in callers && icount_check) {
        clip_qc_args += ' icount '
    }

    if ('paraclu' in callers) {
        clip_qc_args += ' paraclu '
    }

    if ('pureclip' in callers) {
        clip_qc_args += ' pureclip '
    }

    if ('piranha' in callers) {
        clip_qc_args += ' piranha '
    }

    """
    clip_qc.py $clip_qc_args
    """
}

/*
 * STEP 9 - MultiQC
 */
process multiqc {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from ch_fastqc_pretrim_mqc.collect().ifEmpty([])
    file ('cutadapt/*') from ch_cutadapt_mqc.collect().ifEmpty([])
    file ('premap/*') from ch_premap_mqc.collect().ifEmpty([])
    file ('mapped/*') from ch_align_mqc.collect().ifEmpty([])
    path ('preseq/*') from ch_preseq_mqc.collect().ifEmpty([])
    path ('rseqc/*') from ch_rseqc_mqc.collect().ifEmpty([])
    file ('clipqc/*') from ch_clipqc_mqc.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}

/*
 * STEP 10 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/clipseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/clipseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/clipseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/clipseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def max_multiqc_email_size = params.max_multiqc_email_size as nextflow.util.MemoryUnit
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/clipseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/clipseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/clipseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/clipseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

// Check file extension - from nf-core/rnaseq
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
