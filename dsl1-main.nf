#!/usr/bin/env nextflow

//params.input = '/Users/westc/nextflow/dev/data/metadata-char-edit.csv'
//params.star_index = '/Users/westc/nextflow/dev/data/chr20/star_index_2_6'
//ch_star = Channel.value(params.star_index)
// params.bt2_index = ["/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.1.bt2",
// "/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.2.bt2",
// "/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.3.bt2",
// "/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.4.bt2",
// "/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.rev.1.bt2",
// "/Users/westc/nextflow/dev/data/chr20/small_rna_bowtie_ind/small_rna_bowtie_ind.rev.2.bt2"]
// ch_bt2_index = Channel.value(params.bt2_index)
params.fai = '/Users/westc/nextflow/dev/data/chr20/chr20.fa.fai'
ch_fai = Channel.value(params.fai)

opts_cutadapt = params['cutadapt']
opts_bowtie2 = params['bowtie2_align']
opts_star = params['star_align_reads']
opts_umi = params['umi_tools']
opts_crosslinks = params['get_crosslinks']
opts_bt2 = params['bowtie2_align']
opts_fastqc = params['fastqc']

////////////////////////////////////////////////////////////
/* --        iGenomes and smRNA logical actions        -- */
////////////////////////////////////////////////////////////
// Set params and channels based on presence/absence of genome and smRNA associated parameters

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// If igenomes is specified and no smRNA bt2 index is given
if (params.genome && !params.bt2_index && params.genomes.containsKey(params.genome) && params.smRNA.containsKey(params.genome)) {
  params.bt2_index = params.smRNA[ params.genome ].bt2_ind
}


//////////////////////////////////
/* --        iGenomes        -- */
//////////////////////////////////


// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
// params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
// params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
// params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
// params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
ch_star = Channel.value(params.star_index)


///////////////////////////////////////
/* --        Stored smRNA         -- */
///////////////////////////////////////
//params.bt2_index = params.smRNA_species ? params.smRNA[ params.smRNA_species ].bt2_ind ?: false : false
ch_bt2_index = Channel.value(params.bt2_index)


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
ch_fastq_metadata_pre = Channel
                        .fromPath(params.input)
                        .splitCsv(header:true)
                        .map{ row -> processRow(row) }
                        .into{ ch_fastq_metadata;
                               ch_fastqc_pretrim }
                        //.set{ metadata }

//metadata.view()
//def meta = [:]


/*
 * STEP 1.5: FastQC
 */

process fastqc {
    publishDir "${params.outdir}/${opts_fastqc.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts_fastqc.publish_results == "none") null
                      else filename }
    
    container 'biocontainers/fastqc:v0.11.9_cv6'

    input:
        //val opts
        tuple val(meta), path(reads) from ch_fastqc_pretrim

    output:
        path "*.zip", emit: report

    script:
    args = ""
        if(opts_fastqc.args && opts_fastqc.args != '') {
            ext_args = opts_fastqc.args
            args += ext_args.trim()
        }

    prefix = opts_fastqc.suffix ? "${meta.sample_id}${opts_fastqc.suffix}" : "${meta.sample_id}"

    fastqc_command = "fastqc ${args} --threads ${task.cpus} $reads"
    if (params.verbose){
        println ("[MODULE] fastqc command: " + fastqc_command)
    }
    
    //SHELL
    readList = reads.collect{it.toString()}
    if(readList.size > 1){
            """
            ${fastqc_command}
            mv ${reads[0].simpleName}*.zip ${prefix}_fastqc.zip
            mv ${reads[1].simpleName}*.zip ${prefix}_fastqc.zip
            """
    }
    else {
            """
            ${fastqc_command}
            mv ${reads[0].simpleName}*.zip ${prefix}_fastqc.zip
            """
    }
}

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
        tuple val(meta), path(reads) from ch_fastq_metadata

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

process bowtie2_align {
    publishDir "${params.outdir}/${opts_bt2.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts_bt2.publish_results == "none") null
                      else filename }
    
    container 'luslab/nf-modules-bowtie2:latest'

    input:
        //val opts
        tuple val(meta), path(reads) from ch_cut_fastq
        path index from ch_bt2_index

    output:
        tuple val(meta), path("*.sam"), optional: true, emit: sam
        tuple val(meta), path("*.bam"), path("*.bai"), optional: true, emit: bam
        tuple val(meta), path("${prefix}${opts_bt2.unmapped_suffix}.1.fastq.gz"), path("${prefix}${opts_bt2.unmapped_suffix}.2.fastq.gz"), optional: true, emit: unmappedFastqPaired
        tuple val(meta), path("${prefix}${opts_bt2.unmapped_suffix}.fastq.gz") into ch_from_bt2 //, optional: true, emit: unmappedFastqSingle
        path "*stats.txt", emit: report

    script:
        args = "-p ${task.cpus} --no-unal"
        files = ''

        if(opts_bt2.args && opts_bt2.args != '') {
            ext_args = opts_bt2.args
            args += ' ' + ext_args.trim()
        }

        readList = reads.collect{it.toString()}
        if(readList.size > 1){
            files = '-1 ' + reads[0] + ' -2 ' + reads[1]
        }
        else {
            files = '-U ' + reads[0]
        }

        prefix = opts_bt2.suffix ? "${meta.sample_id}${opts_bt2.suffix}" : "${meta.sample_id}"

        // If clause for creating unmapped filename if requested
        if(opts_bt2.unmapped_suffix && opts_bt2.unmapped_suffix != '') {
            if(readList.size > 1){
                args += ' --un-conc-gz ' + "${prefix}${opts_bt2.unmapped_suffix}" + '.fastq.gz'
            }
            else {
                args += ' --un-gz ' + "${prefix}${opts_bt2.unmapped_suffix}" + '.fastq.gz'
            }
        }

        // command = "bowtie2 -x ${index[0].simpleName} $args $files 2>bowtie2_stats.txt > ${prefix}.sam"

        sort_command = "samtools sort -@ ${task.cpus} /dev/stdin > ${prefix}.bam"
        index_command = "samtools index -@ ${task.cpus} ${prefix}.bam"

        if ( opts_bt2.output_sam && opts_bt2.output_sam == true ) {
            command = "bowtie2 -x ${index[0].simpleName} $args $files 2>bowtie2_stats.txt > ${prefix}.sam"
        }
        else {
            command = "bowtie2 -x ${index[0].simpleName} $args $files 2>bowtie2_stats.txt | $sort_command && $index_command"
        }

        if (params.verbose){
            println ("[MODULE] bowtie2 command: " + command)
        }

        """
        $command
        cat bowtie2_stats.txt
        """
}


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
      tuple val(meta), path(reads) from ch_from_bt2//from ch_cut_fastq
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

