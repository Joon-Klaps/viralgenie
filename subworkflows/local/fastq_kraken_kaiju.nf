// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { KRAKEN2_KRAKEN2           } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN           } from '../../../modules/nf-core/bracken/bracken/main'
include { KAIJU_KAIJU               } from '../../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE         } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { UNTAR as UNTAR_KRAKEN2_DB } from '../../../modules/nf-core/kaiju/kaiju/main'
include { UNTAR as UNTAR_BRACKEN_DB } from '../../../modules/nf-core/kaiju/kaiju/main'
include { UNTAR as UNTAR_KAIJU_DB   } from '../../../modules/nf-core/kaiju/kaiju/main'

workflow FASTQ_KRAKEN_KAIJU {

    take:
    reads // channel: [ val(meta), [ fastq ] ]
    kraken2_db
    bracken_db
    kaiju_db

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    ch_versions = Channel.empty()
    if {!params.skip_kraken2}{
            if (kraken2_db.endsWith('.tar.gz')) {
                UNTAR_KRAKEN2_DB (
                    [ [:], kraken2_db ]
                )
                ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_KRAKEN2_DB.out.versions)
            } else {
                        ch_kraken2_db = Channel.value(file(kraken2_db))
            }

        KRAKEN2_KRAKEN2 ( reads, ch_kraken2_db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )

        if (!params.skip_bracken){
            if (bracken_db.endsWith('.tar.gz')) {
                UNTAR_BRACKEN_DB (
                    [ [:], bracken_db ]
                )
                ch_bracken_db = UNTAR_BRACKEN_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_BRACKEN_DB.out.versions)
            } else {
                        ch_bracken_db = Channel.value(file(bracken_db))
            }
            BRACKEN_BRACKEN (KRAKEN2_KRAKEN2.out.report, ch_bracken_db )
            ch_multiqc_files = ch_multiqc_files.mix( BRACKEN_BRACKEN.out.report )
            ch_versions      = ch_versions.mix( BRACKEN_BRACKEN.out.versions.first() )
        }
    }

    if {!params.skip_kaiju}{
        if (kaiju_db.endsWith('.tar.gz')) {
                UNTAR_KAIJU_DB (
                    [ [:], kaiju_db ]
                )
                ch_kaiju_db = UNTAR_KAIJU_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_KAIJU_DB.out.versions)
            } else {
                        ch_kaiju_db = Channel.value(file(kaiju_db))
            }
        KAIJU_KAIJU(reads, ch_kaiju_db)
        ch_versions            = ch_versions.mix( KAIJU_KAIJU.out.versions.first() )

        KAIJU_KAIJU2TABLE ( KAIJU_KAIJU.out.results, ch_input_for_kaiju.db, params.kaiju_taxon_rank)
        ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE.out.summary )

    }

    emit:
    mqc      = ch_multiqc_files   // channel: [ val(meta), multiqc_file ]
    versions = ch_versions        // channel: [ versions.yml ]
}
