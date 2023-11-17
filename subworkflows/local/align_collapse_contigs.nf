include { MAFFT       } from '../../modules/nf-core/mafft/main'
include { MUSCLE      } from '../../modules/nf-core/muscle/main'
include { CAT_CAT     } from '../../modules/nf-core/cat/cat/main'
include { EMBOSS_CONS } from '../../modules/nf-core/emboss/cons/main'

workflow ALIGN_COLLAPSE_CONTIGS {

    take:
    ch_references_members
    aligner

    main:

    if (aligner == "mafft") {
        ch_references      = ch_references_members.map{ it -> [it[0], it[1]] }
        ch_members_no_meta = ch_references_members.map{ it -> it[2] }

        MAFFT ( ch_references, ch_members_no_meta )

        ch_align    = MAFFT.out.fas
        ch_versions = MAFFT.out.versions.first()
    }

    if (aligner == "muscle") {
        ch_sequences = ch_references_members.map{ it -> [it[0], [it[1], it[2]]] }

        CAT_CAT(ch_sequences)

        MUSCLE( CAT_CAT.out.file_out)

        ch_versions = CAT_CAT.out.versions.first()
        ch_align    = MUSLCE.out.aligned_fasta
        ch_versions = ch_versions.mix(MUSCLE.out.versions.first())
    }

    EMBOSS_CONS ( ch_align )
    ch_versions = ch_versions.mix(EMBOSS_CONS.out.versions.first())

    emit:
    consensus      = EMBOSS_CONS.out.consensus // channel: [ val(meta), [ fasta ] ]
    aligned_fasta  = ch_align                  // channel: [ val(meta), [ fasta ] ]
    versions       = ch_versions               // channel: [ versions.yml ]
}

