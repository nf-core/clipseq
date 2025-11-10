include { SAMTOOLS_MERGE                                            } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_GROUP_INDEX                    } from '../../modules/nf-core/samtools/index/main'
include { GET_CROSSLINKS as CALC_CROSSLINKS_INDIVIDUAL              } from '../../modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_CROSSLINKS_INDIVIDUAL_HASGROUP     } from '../../modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_CROSSLINKS_GROUP_HASGROUP          } from '../../modules/local/get_crosslinks'
include { BEDGRAPH_STRAND_SPLIT as STRAND_SPLIT_INDIVIDUAL          } from '../../modules/local/bedgraph_strand_split/main'
include { BEDGRAPH_STRAND_SPLIT as STRAND_SPLIT_INDIVIDUAL_HASGROUP } from '../../modules/local/bedgraph_strand_split/main'
include { BEDGRAPH_STRAND_SPLIT as STRAND_SPLIT_GROUP_HASGROUP      } from '../../modules/local/bedgraph_strand_split/main'
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_POS_INDIVIDUAL            } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_NEG_INDIVIDUAL            } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_POS_INDIVIDUAL_HASGROUP   } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_NEG_INDIVIDUAL_HASGROUP   } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_POS_GROUP_HASGROUP        } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_NEG_GROUP_HASGROUP        } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"

workflow RESOLVE_GROUPS_AND_CROSSLINKS {
    take:
    ch_bam         // channel: [ val(meta), [ bam ] ]
    ch_bai         // channel: [ val(meta), [ bai ] ]
    ch_fasta       // channel: [ val(meta), fasta ]
    ch_fasta_fai   // channel: [ val(meta), fasta.fai ]

    main:
    ch_versions = Channel.empty()

    ch_bam.branch {
        hasGroup: it[0].group  // Branch condition for samples with a group
        noGroup: it[0].group == ''  // Branch condition for samples without a group
    }.set { ch_bam_branches }

    ch_bai.branch {
        hasGroup: it[0].group  // Branch condition for samples with a group
        noGroup: it[0].group == ''  // Branch condition for samples without a group
    }.set { ch_bai_branches }

    ch_bam_branches.hasGroup
        .map { item ->
            def meta = item[0]
            def bam = item[1]
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

            return [newMeta, bam]
        }
        .set { ch_grouped_bam }

    SAMTOOLS_MERGE (
        ch_grouped_bam,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    ch_grouped_bam = SAMTOOLS_MERGE.out.bam

    ch_grouped_bam.concat(ch_bam_branches.noGroup).set { ch_peakcalling_bam }

    SAMTOOLS_GROUP_INDEX (
        ch_peakcalling_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_GROUP_INDEX.out.versions)
    ch_peakcalling_bai = SAMTOOLS_GROUP_INDEX.out.bai


    //
    // MODULE: Run crosslink calculation for samples without a group
    //

    CALC_CROSSLINKS_INDIVIDUAL (
        ch_bam_branches.noGroup.join(ch_bai_branches.noGroup),
        ch_fasta_fai
    )
    ch_versions                        = ch_versions.mix(CALC_CROSSLINKS_INDIVIDUAL.out.versions)
    ch_crosslink_INDIVIDUAL_bed        = CALC_CROSSLINKS_INDIVIDUAL.out.bed
    ch_crosslink_INDIVIDUAL_bedgraph   = CALC_CROSSLINKS_INDIVIDUAL.out.bedgraph

    //
    // MODULE: Split bedgraph by strand for samples without a group
    //

    STRAND_SPLIT_INDIVIDUAL (
        ch_crosslink_INDIVIDUAL_bedgraph
    )
    ch_versions = ch_versions.mix(STRAND_SPLIT_INDIVIDUAL.out.versions)

    //
    // MODULE: Convert strand-specific bedgraphs to bigwig for samples without a group
    //

    BIGWIG_POS_INDIVIDUAL (
        STRAND_SPLIT_INDIVIDUAL.out.pos_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_POS_INDIVIDUAL.out.versions)

    BIGWIG_NEG_INDIVIDUAL (
        STRAND_SPLIT_INDIVIDUAL.out.neg_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_NEG_INDIVIDUAL.out.versions)

    //
    // MODULE: Run crosslink calculation for samples WITH A GROUP individually
    //

    CALC_CROSSLINKS_INDIVIDUAL_HASGROUP (
        ch_bam_branches.hasGroup.join(ch_bai_branches.hasGroup),
        ch_fasta_fai
    )
    ch_versions                                 = ch_versions.mix(CALC_CROSSLINKS_INDIVIDUAL_HASGROUP.out.versions)
    ch_crosslink_INDIVIDUAL_HASGROUP_bed        = CALC_CROSSLINKS_INDIVIDUAL_HASGROUP.out.bed
    ch_crosslink_INDIVIDUAL_HASGROUP_bedgraph   = CALC_CROSSLINKS_INDIVIDUAL_HASGROUP.out.bedgraph

    //
    // MODULE: Split bedgraph by strand for samples WITH A GROUP individually
    //

    STRAND_SPLIT_INDIVIDUAL_HASGROUP (
        ch_crosslink_INDIVIDUAL_HASGROUP_bedgraph
    )
    ch_versions = ch_versions.mix(STRAND_SPLIT_INDIVIDUAL_HASGROUP.out.versions)

    //
    // MODULE: Convert strand-specific bedgraphs to bigwig for samples WITH A GROUP individually
    //

    BIGWIG_POS_INDIVIDUAL_HASGROUP (
        STRAND_SPLIT_INDIVIDUAL_HASGROUP.out.pos_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_POS_INDIVIDUAL_HASGROUP.out.versions)

    BIGWIG_NEG_INDIVIDUAL_HASGROUP (
        STRAND_SPLIT_INDIVIDUAL_HASGROUP.out.neg_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_NEG_INDIVIDUAL_HASGROUP.out.versions)

    //
    // MODULE: Run crosslink calculation for samples WITH A GROUP AS A GROUP
    //
    
    CALC_CROSSLINKS_GROUP_HASGROUP (
        ch_grouped_bam.join(ch_peakcalling_bai),
        ch_fasta_fai
    )
    ch_versions                            = ch_versions.mix(CALC_CROSSLINKS_GROUP_HASGROUP.out.versions)
    ch_crosslink_GROUP_HASGROUP_bed        = CALC_CROSSLINKS_GROUP_HASGROUP.out.bed
    ch_crosslink_GROUP_HASGROUP_bedgraph   = CALC_CROSSLINKS_GROUP_HASGROUP.out.bedgraph

    //
    // MODULE: Split bedgraph by strand for samples WITH A GROUP AS A GROUP
    //

    STRAND_SPLIT_GROUP_HASGROUP (
        ch_crosslink_GROUP_HASGROUP_bedgraph
    )
    ch_versions = ch_versions.mix(STRAND_SPLIT_GROUP_HASGROUP.out.versions)

    //
    // MODULE: Convert strand-specific bedgraphs to bigwig for samples WITH A GROUP AS A GROUP
    //

    BIGWIG_POS_GROUP_HASGROUP (
        STRAND_SPLIT_GROUP_HASGROUP.out.pos_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_POS_GROUP_HASGROUP.out.versions)

    BIGWIG_NEG_GROUP_HASGROUP (
        STRAND_SPLIT_GROUP_HASGROUP.out.neg_bedgraph, 
        ch_fasta_fai.map{ it[1] }
    )
    ch_versions = ch_versions.mix(BIGWIG_NEG_GROUP_HASGROUP.out.versions)

    //
    // Combine crosslinking results for moving forwards
    //

    ch_crosslink_bed = ch_crosslink_INDIVIDUAL_bed.mix(ch_crosslink_GROUP_HASGROUP_bed)
    ch_crosslink_combined_pos_bigwig          = BIGWIG_POS_INDIVIDUAL.out.bigwig.mix(BIGWIG_POS_GROUP_HASGROUP.out.bigwig)
    ch_crosslink_combined_neg_bigwig          = BIGWIG_NEG_INDIVIDUAL.out.bigwig.mix(BIGWIG_NEG_GROUP_HASGROUP.out.bigwig)

    emit:
    versions                               = ch_versions
    crosslink_group_resolved               = ch_crosslink_bed                          // channel: [ val(meta), [ bed ] ] 
    crosslink_INDIVIDUAL                   = ch_crosslink_INDIVIDUAL_bed               // channel: [ val(meta), [ bed ] ] 
    crosslink_INDIVIDUAL_HASGROUP          = ch_crosslink_INDIVIDUAL_HASGROUP_bed      // channel: [ val(meta), [ bed ] ] 
    crosslink_GROUP_HASGROUP               = ch_crosslink_GROUP_HASGROUP_bed           // channel: [ val(meta), [ bed ] ]
    
    // Combined BigWigs
    crosslink_combined_pos_bigwig          = ch_crosslink_combined_pos_bigwig          // channel: [ val(meta), [ bigwig ] ]
    crosslink_combined_neg_bigwig          = ch_crosslink_combined_neg_bigwig          // channel: [ val(meta), [ bigwig ] ]

    // Strand-specific BigWig outputs
    crosslink_INDIVIDUAL_pos_bigwig        = BIGWIG_POS_INDIVIDUAL.out.bigwig          // channel: [ val(meta), [ bigwig ] ]
    crosslink_INDIVIDUAL_neg_bigwig        = BIGWIG_NEG_INDIVIDUAL.out.bigwig          // channel: [ val(meta), [ bigwig ] ]
    crosslink_INDIVIDUAL_HASGROUP_pos_bigwig = BIGWIG_POS_INDIVIDUAL_HASGROUP.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    crosslink_INDIVIDUAL_HASGROUP_neg_bigwig = BIGWIG_NEG_INDIVIDUAL_HASGROUP.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    crosslink_GROUP_HASGROUP_pos_bigwig    = BIGWIG_POS_GROUP_HASGROUP.out.bigwig      // channel: [ val(meta), [ bigwig ] ]
    crosslink_GROUP_HASGROUP_neg_bigwig    = BIGWIG_NEG_GROUP_HASGROUP.out.bigwig      // channel: [ val(meta), [ bigwig ] ]
    
    genome_peakcalling_bam                 = ch_peakcalling_bam                        // channel: [ val(meta), [ bam ] ]
    genome_peakcalling_bai                 = ch_peakcalling_bai                        // channel: [ val(meta), [ bai ] ]
}