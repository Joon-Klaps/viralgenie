include { UNTAR     } from '../../modules/nf-core/untar/main'
include { GUNZIP    } from '../../modules/nf-core/gunzip/main'

workflow UNPACK_DB  {

    take:
    db_in  // channel [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()

    db_in
    .branch { meta, db ->
        tar: db.endsWith('.tar.gz') || db.endsWith('.tgz')
        gzip: db.endsWith('.gz')
        dir: file(db).isDirectory() || db.endsWith('.fna') || db.endsWith('.fa')
        other: true
    }
    .set{db}

    if (db.other) {
        db.other.map{ meta, db ->
            Nextflow.error("Database: ${db} not recognized as '.tar.gz', '.tgz' or '.gz' file. \n\n please download and unpack manually and try again specifying the directory.")
        }
    }

    ch_untar = UNTAR(db.tar).untar
    ch_versions   = ch_versions.mix(UNTAR.out.versions)

    ch_gunzip = GUNZIP(db.gzip).gunzip
    ch_versions   = ch_versions.mix(GUNZIP.out.versions)

    ch_db = Channel.empty()
    ch_db = ch_db.mix(db.dir, ch_untar, ch_gunzip)

    emit:
    db       = ch_db            // channel: [ db ]
    versions = ch_versions      // channel: [ versions.yml ]
}

