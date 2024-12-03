
include { filterContigs                                        } from '../../subworkflows/local/utils_nfcore_viralgenie_pipeline'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_SINGLETON } from '../../modules/local/rename_fasta_header'

workflow SINGLETON_FILTERING {

    take:
    fasta               // channel: [ val(meta), [ fasta ] ]
    min_contig_size     // int
    max_n_perc        // int

    main:
    ch_versions = Channel.empty()

    if ( !params.skip_singleton_filtering) {
        filtered = filterContigs ( fasta, min_contig_size, max_n_perc)
        contig = filtered.pass
    }
    // Rename to avoid errors downstream
    RENAME_FASTA_HEADER_SINGLETON(
        contig,
        []
        )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_SINGLETON.out.versions)


    emit:
    filtered     = RENAME_FASTA_HEADER_SINGLETON.out.fasta  // channel: [ val(meta), [ fasta ] ]
    versions     = ch_versions                              // channel: [ versions.yml ]
}

