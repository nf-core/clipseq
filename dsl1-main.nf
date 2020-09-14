#!/usr/bin/env nextflow

//params.input = '/Users/westc/nextflow/dev/data/metadata-char-edit.csv'
params.star_index = '/Users/westc/nextflow/dev/data/chr20/star_index'
ch_star = Channel.fromPath(params.star_index)

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

opts = params['cutadapt']
opts.args = "-j 8 -a AGATCGGAAGAGC -m 12"

process cutadapt {
    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy", 
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
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
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

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
 * STEP 3: Pre-map to rRNA (bowtie, skip for now)
 */



/*
 * STEP 4: Align to Genome
 */
def opts = params['star_align_reads']
    opts.args = "--outFilterMultimapNmax 1 \
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
                --outSAMtype BAM SortedByCoordinate"

process star_align_reads {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    // container 'quay.io/biocontainers/star:2.7.5b--0'
    container 'luslab/luslab-nf-star:latest'

    input:
      // val opts
      tuple val(meta), path(reads) from ch_cut_fastq
      path star_index from ch_star

    output:
      tuple val(meta), path("*.sam"), optional: true, emit: sam_files
      tuple val(meta), path("*.bam"), path("*.bai"), optional: true, emit: bam_files
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
    if ( opts.args ) {
      if ( opts.args =~ /(--solo)/ ) {
        exit 1, "Error: This module does not support STARsolo (--solo* options). For processing of single-cell RNA-seq data with STAR please use a dedicated module. Exit."
      }
      if ( opts.args =~ /(--runMode)/ ) {
        exit 1, "Error: --runMode is automatically set to 'alignReads'. You do not need to provide it manually. Exit."
      }
      if ( opts.args =~ /(--parametersFiles)/ ) {
        exit 1, "Error: Parameter files (--parametersFiles option) are not supported in this module. Please provide all options not covered by input channels and module parameters via the params.modules['star_align_reads'].args parameter. Exit."
      }
      ext_args = opts.args
      args += ext_args.trim() + " "
    }

    prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

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
    if ( opts.sjdbGTFfile ) {
      args += "--sjdbGTFfile ${opts.sjdbGTFfile} "
    }
    if ( opts.sjdbFileChrStartEnd ) {
      args += "--sjdbFileChrStartEnd ${opts.sjdbFileChrStartEnd} "
    }
    if ( opts.varVCFfile ) {
      args += "--varVCFfile ${opts.varVCFfile} "
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


/*
 * STEP 5: PCR Duplicate Removal
 */


/*
 * STEP 6: Get Crosslinks
 */