include { MAFFT       } from '../../../modules/nf-core/mafft/main'
include { MUSCLE      } from '../../../modules/nf-core/muscle/main'
include { CAT_CAT     } from '../../../modules/nf-core/cat/cat/main'
include { EMBOSS_CONS } from '../../../modules/nf-core/emboss/cons/main'

workflow ALIGN_COLLAPSE_CONTIGS {

    take:
    ch_members
    ch_references
    aligner

    main:

    if (aligner == "mafft") {
        ch_align = MAFFT ( ch_members, ch_references ).fas
        ch_versions = MAFFT.out.versions.first()
    }

    if (aligner == "muscle") {
        ch_members_references_joined = ch_members.join(ch_references, remainder: true)
        ch_sequences = ch_members_references_joined.map{ it -> [it[0], [it[1], it[2]]] }
        CAT_CAT(ch_sequences)
        ch_align = MUSCLE( CAT_CAT.out.file_out).aligned_fasta
        ch_versions = MUSCLE.out.versions.first()
    }

    EMBOSS_CONS ( ch_align )
    ch_versions = ch_versions.mix(EMBOSS_CONS.out.versions.first())

    emit:
    consensus      = EMBOSS_CONS.out.consensus // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

