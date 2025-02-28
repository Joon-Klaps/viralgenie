include { UNTAR as UNTAR_DB      } from '../../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_DB    } from '../../modules/nf-core/gunzip/main'
include { XZ_DECOMPRESS as XZ_DB } from '../../modules/nf-core/xz/decompress/main'

workflow UNPACK_DB  {

    take:
    db_in  // channel [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()
    ch_db       = Channel.empty()

    db_in
    .branch { meta, dbs ->
        tar: dbs.name.endsWith('.tar.gz') || dbs.name.endsWith('.tgz') || dbs.name.endsWith('.tar')
        gzip: dbs.name.endsWith('.gz')
        xz: dbs.name.endsWith('.xz')
        other: true
    }
    .set{db}

    ch_db       = ch_db.mix(db.other)

    ch_db       = ch_db.mix(UNTAR_DB(db.tar).untar)
    ch_versions = ch_versions.mix(UNTAR_DB.out.versions.first())

    ch_db       = ch_db.mix(GUNZIP_DB(db.gzip).gunzip)
    ch_versions = ch_versions.mix(GUNZIP_DB.out.versions.first())

    ch_db       = ch_db.mix(XZ_DB(db.xz).file)
    ch_versions = ch_versions.mix(XZ_DB.out.versions.first())

    emit:
    db       = ch_db            // channel: [ db ]
    versions = ch_versions      // channel: [ versions.yml ]
}

