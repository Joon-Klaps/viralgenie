include { CAT_CAT as CAT_CLUSTER                                      } from '../../modules/nf-core/cat/cat/main'
include { MINIMAP2_INDEX as MINIMAP2_CONTIG_INDEX                     } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_CONTIG_ALIGN                     } from '../../modules/nf-core/minimap2/align/main'
include { IVAR_CONSENSUS as IVAR_CONTIG_CONSENSUS                     } from '../../modules/nf-core/ivar/consensus/main'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CONTIG_CONSENSUS } from '../../modules/local/rename_fasta_header'
include { LOWCOV_TO_REFERENCE                                     } from '../../modules/local/lowcov_to_reference/main'

workflow ALIGN_COLLAPSE_CONTIGS {

    take:
    ch_references_members

    main:
    ch_versions = Channel.empty()


    ch_sequences = ch_references_members.map{ meta, references, members -> [meta, [references, members]] }

    CAT_CLUSTER(ch_sequences)
    ch_versions = ch_versions.mix(CAT_CLUSTER.out.versions.first())
    // TODO:
    // if external_reference is false, we need to include the reference when mapping towards the reference (i.e. the reference is also a member)
    // if external_reference is true, we need to exclude the reference when mapping towards the reference
    // We align contigs to the reference using minimap2
    // Call consensus using IVAR_consensus with low treshholds (eg. 1) (needs only 1 coverage)
    // If there are ambigous bases in the consensus due to low coverage we populate it with the reference sequence.

    ch_references = ch_references_members.map{ meta, references, members -> [meta, references] }

    MINIMAP2_CONTIG_INDEX(ch_references)
    ch_versions = ch_versions.mix(MINIMAP2_CONTIG_INDEX.out.versions.first())

    MINIMAP2_CONTIG_INDEX
        .out
        .index
        .join( ch_references_members, by: [0] )
        .join( CAT_CLUSTER.out.file_out, by: [0] )
        .branch{ meta, index, references, members, comb ->
            external: meta.external_reference
                return [meta, index, members]
            internal: true
                return [meta, index, comb]
        }
        .set{ ch_splitup }

    ch_index_contigs = ch_splitup.external.mix(ch_splitup.internal)

    ch_index = ch_index_contigs.map{ meta, index, contigs -> [meta, index] }
    ch_contigs = ch_index_contigs.map{ meta, index, contigs -> [meta, contigs] }

    MINIMAP2_CONTIG_ALIGN(ch_contigs, ch_index, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_CONTIG_ALIGN.out.versions.first())

    ch_references_members
        .join( MINIMAP2_CONTIG_ALIGN.out.bam, by: [0] )
        .map{ meta, references, members, bam -> [meta, references, bam] }
        .set{ ch_references_bam }

    ivar_bam   = ch_references_bam.map{ meta, references, bam -> [meta, bam] }
    ivar_fasta = ch_references_bam.map{ meta, references, bam -> [references] }
    IVAR_CONTIG_CONSENSUS(
        ivar_bam,
        ivar_fasta,
        true
    )
    ch_versions= ch_versions.mix(IVAR_CONTIG_CONSENSUS.out.versions.first())

    RENAME_FASTA_HEADER_CONTIG_CONSENSUS( IVAR_CONTIG_CONSENSUS.out.fasta, "consensus" )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_CONTIG_CONSENSUS.out.versions.first())

    // If external, there possibly regions that require patching
    RENAME_FASTA_HEADER_CONTIG_CONSENSUS.out.fasta
        .branch{ meta, fasta ->
            external: meta.external_reference
            internal: true
        }
        .set{ ch_consensus }

    // Combine input for custom annotation script.
    ch_references_bam
        .join( ch_consensus.external, by: [0] )
        .join( IVAR_CONTIG_CONSENSUS.out.mpileup, by:[0])
        .map{ meta, references, bam, consensus, mpileup -> [meta, references, consensus, mpileup] }
        .set{ ch_ref_cons_mpileup }


    aligned_txt = Channel.empty()
    // Custom script that replaces region in consensus with orignally 0 coverage with regions from the reference.
    if (!params.skip_hybrid_consensus) {
        LOWCOV_TO_REFERENCE( ch_ref_cons_mpileup)
        ch_versions = ch_versions.mix(LOWCOV_TO_REFERENCE.out.versions.first())
        aligned_txt = LOWCOV_TO_REFERENCE.out.txt

        LOWCOV_TO_REFERENCE
            .out
            .sequence
            .mix( ch_consensus.internal )
            .set{ consensus_patched }
    } else {
        consensus_patched = RENAME_FASTA_HEADER_CONTIG_CONSENSUS.out.fasta
    }

    emit:
    consensus       = consensus_patched                 // channel: [ val(meta), [ fasta ] ]
    aligned_txt     = aligned_txt                       // channel: [ val(meta), [ txt ] ]
    unaligned_fasta = CAT_CLUSTER.out.file_out          // channel: [ val(meta), [ fasta ] ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}

