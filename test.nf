#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Simulate input data
    def inputData = Channel.of(
        [[id:'PHO92_A', group:'PHO92', control:'ctrl1', single_end:true], '/path/to/PHO92_A.genome.dedup.bam'],
        [[id:'PHO92_B', group:'PHO92', control:'ctrl2', single_end:true], '/path/to/PHO92_B.genome.dedup.bam'],
        [[id:'TDP_A', group:'', control:'ctrl2', single_end:true], '/path/to/PHO92_B.genome.dedup.bam'],
        [[id:'PHNRNPX', group:'', control:'ctrl2', single_end:true], '/path/to/PHO92_B.genome.dedup.bam']
        // Add more simulated items if needed
    )
    inputData.branch {
        hasGroup: it[0].group  // Branch condition for samples with a group
        noGroup: it[0].group == ''  // Branch condition for samples without a group
    }.set { ch_branches }  // Capture branching result into ch_branches
    
    ch_branches.noGroup.view { item -> println("No group item: $item") }
    ch_branches.hasGroup.view { item -> println("Has group item: $item") }

    // Process the input data similarly to your main pipeline logic
    ch_branches.hasGroup
        .map { item ->
            println("Mapping item: $item")
            def meta = item[0]
            def bam = item[1]
            println("Processed mapping: Meta=${meta}, BAM=${bam}")
            return [meta.group, meta, bam]
        }
        .groupTuple(by: 0)
        .map { tuple ->
            def group = tuple[0]
            def items = tuple[1]
            def bam = tuple[2]

            def newMeta = [:]
            newMeta.id = group
            newMeta.group = group
            newMeta.control = items[0].control
            newMeta.single_end = true

            println("NewMeta: $newMeta")

            return [newMeta, bam]
        }
        .view { item -> println("Final Processed item: $item") }
}