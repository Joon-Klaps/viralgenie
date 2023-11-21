include { MAFFT       } from '../../modules/nf-core/mafft/main'
include { MUSCLE      } from '../../modules/nf-core/muscle/main'
include { CAT_CAT     } from '../../modules/nf-core/cat/cat/main'
include { EMBOSS_CONS } from '../../modules/nf-core/emboss/cons/main'

workflow ALIGN_COLLAPSE_CONTIGS {

    take:
    ch_references_members
    aligner

    main:

    ch_sequences = ch_references_members.map{ it -> [it[0], [it[1], it[2]]] }

    CAT_CAT(ch_sequences)

    if (aligner == "mafft") {
        ch_references      = ch_references_members.map{ meta,centroid,members -> [meta, centroid] }
        ch_members         = ch_references_members.map{ meta,centroid,members -> [meta, members] }

        MAFFT (
            ch_references,
            [[:],[]],
            ch_members,
            [[:],[]],
            [[:],[]],
            [[:],[]]
            )

        ch_align    = MAFFT.out.fas
        ch_versions = MAFFT.out.versions.first()
    }

    if (aligner == "muscle") {

        MUSCLE( CAT_CAT.out.file_out)

        ch_versions = CAT_CAT.out.versions.first()
        ch_align    = MUSLCE.out.aligned_fasta
        ch_versions = ch_versions.mix(MUSCLE.out.versions.first())
    }

    EMBOSS_CONS ( ch_align )
    ch_versions = ch_versions.mix(EMBOSS_CONS.out.versions.first())

    emit:
    consensus       = EMBOSS_CONS.out.consensus // channel: [ val(meta), [ fasta ] ]
    aligned_fasta   = ch_align                  // channel: [ val(meta), [ fasta ] ]
    unaligned_fasta = CAT_CAT.out.file_out      // channel: [ val(meta), [ fasta ] ]
    versions        = ch_versions               // channel: [ versions.yml ]
}

