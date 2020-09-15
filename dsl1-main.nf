#!/usr/bin/env nextflow

//params.input = '/Users/westc/nextflow/dev/data/metadata-char-edit.csv'
params.star_index = '/Users/westc/nextflow/dev/data/chr20/star_index'
ch_star = Channel.fromPath(params.star_index)
params.fai = '/Users/westc/nextflow/dev/data/chr20/chr20.fa.fai'
ch_fai = Channel.fromPath(params.fai)
params.bt2_index = '/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind'


opts_cutadapt = params['cutadapt']
opts_bowtie2 = params['bowtie2_align']
opts_star = params['star_align_reads']
opts_umi = params['umi_tools']
opts_crosslinks = params['get_crosslinks']

//// Option params ////
//params.opts.cutadapt = params['cutadapt']
//params.opts.cutadapt.args = "-j 8 -a AGATCGGAAGAGC -m 12"
//////////////////////

/////// Process needed for Step 1 /////////
def processRow(LinkedHashMap row) {
    def meta = [:]
    meta.sample_id = row.sample_id

    for (Map.Entry<String, ArrayList<String>> entry : row.entrySet()) {
        String key = entry.getKey();
        String value = entry.getValue();
    
        if(key != "sample_id" && key != "data1" && key != "data2") {
            meta.put(key, value)
        }
    }

    def array = []
    if (row.data2 == null) {
        array = [ meta, [ file(row.data1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.data1, checkIfExists: true), file(row.data2, checkIfExists: true) ] ]
    }
    return array
}
//////////////////////////////////////////

/*
 * STEP 1: FastQ
 */
// setting up metadata structure
ch_fastq_metadata = Channel
                        .fromPath(params.input)
                        .splitCsv(header:true)
                        .map{ row -> processRow(row) }
                        .set{ metadata }

//metadata.view()

//def meta = [:]



/*
 * STEP 2: Read Trimming
 */

//opts = params['cutadapt']
//opts.args = "-j 8 -a AGATCGGAAGAGC -m 12"

process cutadapt {
    publishDir "${params.outdir}/${opts_cutadapt.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (opts_cutadapt.publish_results == "none") null
                    else filename }
    
    container 'quay.io/biocontainers/cutadapt:2.10--py37hf01694f_1'

    input:
        //val opts
        tuple val(meta), path(reads) from metadata

    output:
        tuple val(meta), path("*.fq.gz") into ch_cut_fastq//, emit: fastq
        path "*.log", emit: report

    script:

        args = ""
        if(opts_cutadapt.args && opts_cutadapt.args != '') {
            ext_args = opts_cutadapt.args
            args += ext_args.trim()
        }

        prefix = opts_cutadapt.suffix ? "${meta.sample_id}${opts_cutadapt.suffix}" : "${meta.sample_id}"

        // Construct CL line
        readList = reads.collect{it.toString()}
        if (readList.size > 1){
            cutadapt_command = "cutadapt ${args} -o ${prefix}.1.fq.gz -p ${prefix}.2.fq.gz $reads > ${meta.sample_id}_cutadapt.log"
        } else {
            cutadapt_command = "cutadapt ${args} -o ${prefix}.fq.gz $reads > ${meta.sample_id}_cutadapt.log"
        }
        // Log
        if (params.verbose){
            println ("[MODULE] cutadapt command: " + cutadapt_command)
        }

        //SHELL
        """
        ${cutadapt_command}
        """
}




/*
 * STEP 3: Pre-map to rRNA (bowtie2)
 */



/*
 * STEP 4: Align to Genome
 */

process star_align_reads {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts_star.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts_star.publish_results == "none") null
                      else filename }

    // container 'quay.io/biocontainers/star:2.7.5b--0'
    container 'luslab/luslab-nf-star:latest'

    input:
      // val opts
      tuple val(meta), path(reads) from ch_cut_fastq
      path star_index from ch_star

    output:
      tuple val(meta), path("*.sam"), optional: true, emit: sam_files
      tuple val(meta), path("*.bam"), path("*.bai") into ch_aligned //, optional: true, emit: bam_files
      tuple val(meta), path("*.SJ.out.tab"), optional: true, emit: sj_files
      tuple val(meta), path("*.junction"), optional: true, emit: ch_junctions
      tuple val(meta), path("*.ReadsPerGene.out.tab"),  optional: true, emit: reads_per_gene
      tuple val(meta), path("*.Log.final.out"), emit: final_log_files
      tuple val(meta), path("*.Log.out"), emit: out_log_files
      tuple val(meta), path("*.Log.progress.out"), emit: progress_log_files
      path "*.Log.final.out", emit: report

    script:

    // Add the main arguments
    args = "--runMode alignReads --runDirPerm All_RWX --genomeDir $star_index --readFilesIn $reads "

    // Check and add custom arguments
    if ( opts_star.args ) {
      if ( opts_star.args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( opts_star.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( opts_star.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the params.modules['star_align_reads'].args parameter. Exit."
      }
      ext_args = opts_star.args
      args += ext_args.trim() + " "
    }

    prefix = opts_star.suffix ? "${meta.sample_id}${opts_star.suffix}" : "${meta.sample_id}"

    // Add the number of threads
    args += "--runThreadN $task.cpus "

    // Add output file name prefix
    args += "--outFileNamePrefix ${prefix}. "

    // Add compression parameters 
    test_file_name = "$reads"
    if ( "$test_file_name" =~ /(.gz$)/ ) {
      args += "--readFilesCommand gunzip -c "
    } 
    if ( "$test_file_name" =~ /(.bz2$)/ ) {
      args += "--readFilesCommand bunzip2 -c "
    }

    // Add optional input parameters
    if ( opts_star.sjdbGTFfile ) {
      args += "--sjdbGTFfile ${opts_star.sjdbGTFfile} "
    }
    if ( opts_star.sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${opts_star.sjdbFileChrStartEnd} "
    }
    if ( opts_star.varVCFfile ) {
      args += "--varVCFfile ${opts_star.varVCFfile} "
    }

    // Add memory constraints
    avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000} " : ''
    avail_mem += task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    args += avail_mem

    index_command = "samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam"

    // Construct command line
    map_command = "STAR $args && $index_command"

    // Log
    if (params.verbose) {
        println ("[MODULE] star_align_reads command: " + map_command)
    }

    // Run read mapping with STAR
    """
    ${map_command}
    """
}




// /*
//  * STEP 5: PCR Duplicate Removal
//  */

process umitools_dedup {
    publishDir "${params.outdir}/${opts_umi.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts_umi.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-umi_tools:latest'

    input:
        //val opts
        tuple val(meta), path(bam), path(bai) from ch_aligned
       
    output:
        tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai") into ch_dup_removed//, emit: dedupBam
        path "*.log", emit: report

    script:

    // Init
    prefix = opts_umi.suffix ? "${meta.sample_id}${opts_umi.suffix}" : "${meta.sample_id}"

    args = "--log=${prefix}.log "

    if(opts_umi.args && opts_umi.args != '') {
        ext_args = opts_umi.args
        args += ext_args.trim()
    }

    // Construct CL line
    dedup_command = "umi_tools dedup ${args} -I ${bam[0]} -S ${prefix}.bam --output-stats=${prefix}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    samtools index ${prefix}.bam
    """
}


// /*
//  * STEP 6: Get Crosslinks
//  */

process getcrosslinks {
    publishDir "${params.outdir}/${opts_crosslinks.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts_crosslinks.publish_results == "none") null
                      else filename }

    container 'luslab/nf-modules-get_crosslinks:latest'

    input:
      //val(opts)
      tuple val(meta), path(bam), path(bai) from ch_dup_removed
      path(fai) from ch_fai

    output:
      tuple val(meta), path ("${prefix}.bed.gz"), emit: crosslinkBed

    script:

      prefix = opts_crosslinks.suffix ? "${meta.sample_id}${opts_crosslinks.suffix}" : "${meta.sample_id}"

      //SHELL
      """
      bedtools bamtobed -i ${bam[0]} > dedupe.bed
      bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
      bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
      bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
      cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${prefix}.bed.gz
      """
}

