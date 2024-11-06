//
// Determine metagenomic diversity using Kraken2 and Kaiju
//
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA      } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN                } from '../../modules/nf-core/bracken/bracken/main'
include { KAIJU_KAIJU                    } from '../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE              } from '../../modules/nf-core/kaiju/kaiju2table/main'
include { KAIJU_KAIJU2KRONA              } from '../../modules/nf-core/kaiju/kaiju2krona/main'
include { KRONA_CLEANUP                  } from '../../modules/local/krona_cleanup/main'
include { KRONA_KTIMPORTTEXT             } from '../../modules/nf-core/krona/ktimporttext/main'

workflow FASTQ_KRAKEN_KAIJU {

    take:
    reads            // channel: [ val(meta), [ fastq ] ]
    kraken2_db       // channel: [ path(kraken2_db) ]
    bracken_db       // channel: [ path(bracken_db) ]
    kaiju_db         // channel: [ path(kaiju_db) ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_krona_text           = Channel.empty()
    ch_raw_classifications  = Channel.empty()
    read_classifiers        = params.read_classifiers ? params.read_classifiers.split(',').collect{ it.trim().toLowerCase() } : []


    // Kraken
    if ('kraken2' in read_classifiers){
        KRAKEN2_KRAKEN2 ( reads, kraken2_db, params.kraken2_save_reads, params.kraken2_save_readclassification )
        ch_raw_classifications = ch_raw_classifications.mix(KRAKEN2_KRAKEN2.out.classified_reads_assignment)
        kraken2_report         = KRAKEN2_KRAKEN2.out.report.map{ meta, report -> [meta + [tool: 'kraken2'], report]}
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )

        // Bracken: get more accurate estimates of abundance
        if ('bracken' in read_classifiers){
            BRACKEN_BRACKEN ( kraken2_report, bracken_db )
            ch_versions    = ch_versions.mix( BRACKEN_BRACKEN.out.versions.first() )
            kraken2_report = BRACKEN_BRACKEN.out.reports.map{ meta, report -> [meta + [tool: 'bracken'], report]}
        }


        KRAKENTOOLS_KREPORT2KRONA ( kraken2_report )
        ch_krona_text     = ch_krona_text.mix( KRAKENTOOLS_KREPORT2KRONA.out.txt )
        ch_versions       = ch_versions.mix( KRAKENTOOLS_KREPORT2KRONA.out.versions.first() )
        ch_multiqc_files  = ch_multiqc_files.mix( kraken2_report )
    }

    // Kaiju
    if ('kaiju' in read_classifiers){
        KAIJU_KAIJU(reads, kaiju_db)
        kaiju_report     = KAIJU_KAIJU.out.results.map{ meta, report -> [meta + [tool: 'kaiju'], report]}
        ch_versions      = ch_versions.mix( KAIJU_KAIJU.out.versions.first() )

        KAIJU_KAIJU2TABLE ( kaiju_report, kaiju_db, params.kaiju_taxon_rank)
        ch_multiqc_files = ch_multiqc_files.mix( KAIJU_KAIJU2TABLE.out.summary )
        ch_versions      = ch_versions.mix( KAIJU_KAIJU2TABLE.out.versions )

        KAIJU_KAIJU2KRONA( kaiju_report, kaiju_db )
        ch_krona_text = ch_krona_text.mix( KAIJU_KAIJU2KRONA.out.txt )
        ch_versions   = ch_versions.mix( KAIJU_KAIJU2KRONA.out.versions.first() )

    }

    /*
        Remove taxonomy level annotation from the Krona text files
    */
    KRONA_CLEANUP( ch_krona_text )
    ch_cleaned_krona_text = KRONA_CLEANUP.out.txt
    ch_versions           = ch_versions.mix( KRONA_CLEANUP.out.versions.first() )

    /*
        Convert Krona text files into html Krona visualizations
    */
    ch_krona_text_for_import = ch_cleaned_krona_text
        .map{meta, txt -> [[id: meta.tool], txt]}
        .groupTuple()

    KRONA_KTIMPORTTEXT( ch_krona_text_for_import )
    ch_krona_html = KRONA_KTIMPORTTEXT.out.html
    ch_versions   = ch_versions.mix( KRONA_KTIMPORTTEXT.out.versions.first() )

    emit:
    read_classifications = ch_raw_classifications // channel: [ val(meta), [ classified_reads ] ]
    krona_html           = ch_krona_html          // channel: [ val(meta), [ html ] ]
    mqc                  = ch_multiqc_files       // channel: [ val(meta), multiqc_file ]
    versions             = ch_versions            // channel: [ versions.yml ]
}

