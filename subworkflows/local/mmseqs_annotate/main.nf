include { MMSEQS_CREATEDB as MMSEQS_CREATEANNOTATIONDB } from '../../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_EASYSEARCH                            } from '../../../modules/nf-core/mmseqs/easysearch/main'

workflow MMSEQS_ANNOTATE {
    take:
    ch_genomes // channel: [ val(meta), [ fasta ] ]
    ch_db      // channel: [ val(meta), [ fasta ] ]

    main:

    // create mmseqs annotation db
    MMSEQS_CREATEANNOTATIONDB(ch_db)

    // search the genomes against the annotation db
    MMSEQS_EASYSEARCH(ch_genomes, MMSEQS_CREATEANNOTATIONDB.out.db)

    emit:
    tsv      = MMSEQS_EASYSEARCH.out.tsv // channel: [ val(meta), [ tsv ] ]
}
