include { filterContigs                                        } from '../utils_nfcore_viralmetagenome_pipeline'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_SINGLETON } from '../../../modules/local/rename_fasta_header/main'

workflow SINGLETON_FILTERING {

    take:
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    min_contig_size   // int
    max_n_perc        // int

    main:

    if ( !params.skip_singleton_filtering) {
        ch_filtered = filterContigs ( ch_fasta, min_contig_size, max_n_perc)
        ch_contig   = ch_filtered.pass
    }
    // Rename to avoid errors downstream
    RENAME_FASTA_HEADER_SINGLETON(
        ch_contig,
        []
        )


    emit:
    filtered     = RENAME_FASTA_HEADER_SINGLETON.out.fasta  // channel: [ val(meta), [ fasta ] ]
}
