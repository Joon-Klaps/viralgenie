//
// Determine metagenomic diversity using Kraken2 and Kaiju
//
include { KRAKEN2_KRAKEN2           } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN           } from '../../../modules/nf-core/bracken/bracken/main'
include { KAIJU_KAIJU               } from '../../../modules/nf-core/kaiju/kaiju/main'
include { KAIJU_KAIJU2TABLE         } from '../../../modules/nf-core/kaiju/kaiju2table/main'
include { KAIJU_KAIJU2KRONA         } from '../../../modules/nf-core/kaiju/kaiju2krona/main'
include { KRONA_CLEANUP             } from '../../../modules/local/krona_cleanup/main'
include { KRONA_KTIMPORTTEXT        } from '../../../modules/nf-core/krona/ktimporttext/main'

workflow FASTQ_KRAKEN_KAIJU {
    take:
    ch_reads         // channel: [ val(meta), [ fastq ] ]
    read_classifiers // value ['kraken2','kaiju','bracken']
    ch_kraken2_db    // channel: [ path(kraken2_db) ]
    ch_bracken_db    // channel: [ path(bracken_db) ]
    ch_kaiju_db      // channel: [ path(kaiju_db) ]
    kraken2_save_reads // value: true/false
    kraken2_save_readclassification // value: true/false
    kaiju_taxon_rank // value: phylum, class, order, family, genus, species

    main:
    ch_multiqc_files       = channel.empty()
    ch_krona_text          = channel.empty()
    ch_raw_classifications = channel.empty()


    // Kraken
    if ('kraken2' in read_classifiers || 'bracken' in read_classifiers) {
        KRAKEN2_KRAKEN2(ch_reads, ch_kraken2_db, kraken2_save_reads, kraken2_save_readclassification)
        ch_raw_classifications = ch_raw_classifications.mix(KRAKEN2_KRAKEN2.out.classified_reads_assignment)
        kraken2_report = KRAKEN2_KRAKEN2.out.report.map { meta, report -> [meta + [tool: 'kraken2'], report] }

        // Bracken: get more accurate estimates of abundance, can only run after kraken2
        if ('bracken' in read_classifiers) {
            BRACKEN_BRACKEN(kraken2_report, ch_bracken_db)
            kraken2_report = BRACKEN_BRACKEN.out.txt.map { meta, report -> [meta + [tool: 'bracken'], report] }
        }

        KRAKENTOOLS_KREPORT2KRONA(kraken2_report)
        ch_krona_text    = ch_krona_text.mix(KRAKENTOOLS_KREPORT2KRONA.out.txt)
        ch_multiqc_files = ch_multiqc_files.mix(kraken2_report)
    }

    // Kaiju
    if ('kaiju' in read_classifiers) {
        KAIJU_KAIJU(ch_reads, ch_kaiju_db)
        ch_kaiju_report = KAIJU_KAIJU.out.results.map { meta, report -> [meta + [tool: 'kaiju'], report] }

        KAIJU_KAIJU2TABLE(ch_kaiju_report, ch_kaiju_db, kaiju_taxon_rank)
        ch_multiqc_files = ch_multiqc_files.mix(KAIJU_KAIJU2TABLE.out.summary)

        KAIJU_KAIJU2KRONA(ch_kaiju_report, ch_kaiju_db)
        ch_krona_text = ch_krona_text.mix(KAIJU_KAIJU2KRONA.out.txt)
    }

    /*
        Remove taxonomy level annotation from the Krona text files
    */
    KRONA_CLEANUP(ch_krona_text)
    ch_cleaned_krona_text = KRONA_CLEANUP.out.txt

    /*
        Convert Krona text files into html Krona visualizations
    */
    ch_krona_text_for_import = ch_cleaned_krona_text
        .map { meta, txt -> [[id: meta.tool], txt] }
        .groupTuple()

    KRONA_KTIMPORTTEXT(ch_krona_text_for_import)
    ch_krona_html = KRONA_KTIMPORTTEXT.out.html

    emit:
    read_classifications = ch_raw_classifications // channel: [ val(meta), [ classified_reads ] ]
    krona_html           = ch_krona_html // channel: [ val(meta), [ html ] ]
    mqc                  = ch_multiqc_files // channel: [ val(meta), multiqc_file ]
}
