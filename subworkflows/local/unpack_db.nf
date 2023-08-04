include { UNTAR     } from '../../modules/nf-core/untar/main'
include { GUNZIP    } from '../../modules/nf-core/gunzip/main'

workflow UNPACK_DB  {

    take:
    db  // path: [ db ]

    main:
    ch_versions = Channel.empty()

    if (db.endsWith('.tar.gz')){
        ch_db = UNTAR([ [:], db ]).untar
        ch_versions   = ch_versions.mix(UNTAR.out.versions)
    } else if(db.endsWith('.gz')){
        ch_db = GUNZIP([ [:], db ]).gunzip
        ch_versions   = ch_versions.mix(GUNZIP.out.versions)
    } else if (db){
        ch_db = Channel.value(file(db))
    } else (
        Nextflow.error("Database: ${db} not recognized as '.tar.gz' or '.gz' file. \n\n please unpack and try again specifying the directory.")
    )

    emit:
    db = ch_db                  // channel: [ db ]
    versions = ch_versions      // channel: [ versions.yml ]
}

