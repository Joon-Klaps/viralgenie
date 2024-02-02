
include { MMSEQS_CREATEDB    as MMSEQS_CREATEANNOTATIONDB  } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_CREATEINDEX                               } from '../../modules/nf-core/mmseqs/createindex/main'
include { MMSEQS_EASYSEARCH                                } from '../../modules/nf-core/mmseqs/easysearch/main'

workflow MMSEQS_ANNOTATE {

    take:
    genomes // channel: [ val(meta), [ fasta ] ]
    db      // channel: [ val(meta), [ fasta ] ]
    main:
    ch_versions = Channel.empty()

    // create mmseqs annotation db
    MMSEQS_CREATEANNOTATIONDB ( db )
    ch_versions = ch_versions.mix(MMSEQS_CREATEANNOTATIONDB.out.versions.first())

    // create mmseqs index for faster reuse of the db
    MMSEQS_CREATEINDEX ( MMSEQS_CREATEANNOTATIONDB.out.db )
    ch_versions = ch_versions.mix(MMSEQS_CREATEINDEX.out.versions.first())

    // search the genomes against the annotation db
    MMSEQS_EASYSEARCH ( genomes, MMSEQS_CREATEINDEX.out.db_indexed )
    ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    emit:
    tsv      = MMSEQS_EASYSEARCH.out.tsv // channel: [ val(meta), [ bam ] ]
    versions = ch_versions               // channel: [ versions.yml ]
}

