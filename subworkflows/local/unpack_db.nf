include { UNTAR as UNTAR_DB     } from '../../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_DB    } from '../../modules/nf-core/gunzip/main'

workflow UNPACK_DB  {

    take:
    db_in  // channel [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()

    db_in
    .branch { meta, db ->
        tar: db.name.endsWith('.tar.gz') || db.name.endsWith('.tgz') || db.name.endsWith('.tar')
        gzip: db.name.endsWith('.gz')
        other: true
    }
    .set{db}

    ch_untar = UNTAR_DB(db.tar).untar
    ch_versions   = ch_versions.mix(UNTAR_DB.out.versions.first())

    ch_gunzip = GUNZIP_DB(db.gzip).gunzip
    ch_versions   = ch_versions.mix(GUNZIP_DB.out.versions.first())

    ch_db = Channel.empty()
    ch_db = ch_db.mix(db.other, ch_untar, ch_gunzip)

    emit:
    db       = ch_db            // channel: [ val(meta), [ db ] ]
    versions = ch_versions      // channel: [ versions.yml ]
}

