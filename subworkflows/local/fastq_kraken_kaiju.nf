//
// Determine metagenomic diversity using Kraken2 and Kaiju
//
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN                } from '../../modules/nf-core/bracken/bracken/main'
include { KAIJU_KAIJU                    } from '../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE              } from '../../modules/nf-core/kaiju/kaiju2table/main'

workflow FASTQ_KRAKEN_KAIJU {

    take:
    reads // channel: [ val(meta), [ fastq ] ]
    kraken2_db
    bracken_db
    kaiju_db

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_raw_classifications  = Channel.empty()

    // Kraken
    if (!params.skip_kraken2){
        // decompress kraken2_db if needed
        KRAKEN2_KRAKEN2 ( reads, ch_kraken2_db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_raw_classifications = ch_raw_classifications.mix(KRAKEN2_KRAKEN2.out.classified_reads_assignment)
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )

        // Bracken: get more accurate estimates of abundance
        // TODO TEST
        if (!params.skip_bracken){
            // decompress bracken_db if needed
            BRACKEN_BRACKEN (KRAKEN2_KRAKEN2.out.report, ch_bracken_db )
            ch_versions   = ch_versions.mix( BRACKEN_BRACKEN.out.versions.first() )
        }
    }

    // Kaiju
    if (!params.skip_kaiju){
        // decompress kaiju_db if needed
        KAIJU_KAIJU(reads, ch_kaiju_db)
        ch_versions = ch_versions.mix( KAIJU_KAIJU.out.versions.first() )

        KAIJU_KAIJU2TABLE ( KAIJU_KAIJU.out.results, ch_kaiju_db, params.kaiju_taxon_rank)
        ch_versions = ch_versions.mix( KAIJU_KAIJU2TABLE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE.out.summary )

    }

    emit:
    read_classifications = ch_raw_classifications // channel: [ val(meta), [ classified_reads ] ]
    mqc                  = ch_multiqc_files       // channel: [ val(meta), multiqc_file ]
    versions             = ch_versions            // channel: [ versions.yml ]
}

