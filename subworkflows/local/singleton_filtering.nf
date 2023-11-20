
include { FASTA_CONTIG_FILTERING                               } from './fasta_contig_filtering'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_SINGLETON } from '../../modules/local/rename_fasta_header'

workflow SINGLETON_FILTERING {

    take:
    fasta               // channel: [ val(meta), [ fasta ] ]
    min_contig_size     // int
    max_n_1OOkbp        // int

    main:
    ch_versions = Channel.empty()

    // Rename to avoid errors downstream
    RENAME_FASTA_HEADER_SINGLETON(
        fasta,
        "singleton.contig"
        )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_SINGLETON.out.versions)

    FASTA_CONTIG_FILTERING (
        RENAME_FASTA_HEADER_SINGLETON.out.fasta,
        min_contig_size,
        max_n_1OOkbp
        )
    ch_versions = ch_versions.mix(FASTA_CONTIG_FILTERING.out.versions)

    emit:
    filtered     = FASTA_CONTIG_FILTERING.out.contigs           // channel: [ val(meta), [ fasta ] ]
    renamed      = RENAME_FASTA_HEADER_SINGLETON.out.fasta      // channel: [ val(meta), [ fasta ] ]

    versions     = ch_versions                                  // channel: [ versions.yml ]
}

