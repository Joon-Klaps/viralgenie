include { UNTAR     } from '../../modules/nf-core/untar/main'
include { GUNZIP    } from '../../modules/nf-core/gunzip/main'

workflow UNPACK_DB  {

    take:
    db_in  // channel [ [ db ] ]

    main:
    ch_versions = Channel.empty()

    db_in
    .branch { db ->
        tar: db.name.endsWith('.tar.gz') || db.name.endsWith('.tgz')
            return [[id:db.name], db]
        gzip: db.name.endsWith('.gz')
            return [[id:db.name], db]
        other: true
            return [[id:db.name], db]
    }
    .set{db}

    ch_untar = UNTAR(db.tar).untar
    ch_versions   = ch_versions.mix(UNTAR.out.versions)

    ch_gunzip = GUNZIP(db.gzip).gunzip
    ch_versions   = ch_versions.mix(GUNZIP.out.versions)

    ch_db = Channel.empty()
    ch_db = ch_db.mix(db.other, ch_untar, ch_gunzip)


    emit:
    db       = ch_db            // channel: [ db ]
    versions = ch_versions      // channel: [ versions.yml ]
}

