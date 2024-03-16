
include { filterContigs                                        } from '../../modules/local/functions'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_SINGLETON } from '../../modules/local/rename_fasta_header'

workflow SINGLETON_FILTERING {

    take:
    fasta               // channel: [ val(meta), [ fasta ] ]
    min_contig_size     // int
    max_n_perc        // int

    main:
    ch_versions = Channel.empty()

    contig = filterContigs ( fasta, min_contig_size, max_n_perc )

    // Rename to avoid errors downstream
    RENAME_FASTA_HEADER_SINGLETON(
        contig.pass,
        "singleton.contig"
        )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_SINGLETON.out.versions)


    emit:
    filtered     = contig.pass                              // channel: [ val(meta), [ fasta ] ]
    renamed      = RENAME_FASTA_HEADER_SINGLETON.out.fasta  // channel: [ val(meta), [ fasta ] ]
    versions     = ch_versions                              // channel: [ versions.yml ]
}

