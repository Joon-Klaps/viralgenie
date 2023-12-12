include { KAIJU_KAIJU     as KAIJU_CONTIG      } from '../../modules/nf-core/kaiju/kaiju/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG    } from '../../modules/nf-core/kraken2/kraken2/main'
include { KAIJU_MERGEOUTPUTS                   } from '../../modules/local/kaiju/mergeoutputs/main'

workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs         // channel: [ val(meta), [ fasta ] ]
    ch_kaiju_db        // channel: [ val(meta), [ db ] ]
    ch_kraken2_db      // channel: [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()

    KAIJU_CONTIG ( ch_contigs, ch_kaiju_db)
    kaiju       = KAIJU_CONTIG.out.results
    ch_versions = ch_versions.mix( KAIJU_CONTIG.out.versions.first() )

    KRAKEN2_CONTIG ( ch_contigs, ch_kraken2_db, false, true )
    // combine classified & non classified reads into one channel which we merge later on
    kraken       = KRAKEN2_CONTIG.out.classified_reads_fastq.join(KRAKEN2_CONTIG.out.unclassified_reads_fastq, by: [0], remainder: true)
    ch_versions  = ch_versions.mix( KRAKEN2_CONTIG.out.versions.first() )

    ch_classifications = kaiju.join(kraken) // channel: [ val(meta), [ kaiju ] , [ classified ], [ unclassified ] ]

    KAIJU_MERGEOUTPUTS ( ch_classifications, ch_kaiju_db )
    merged_classifications = KAIJU_MERGEOUTPUTS.out.merged
    ch_versions            = ch_versions.mix( KAIJU_MERGEOUTPUTS.out.versions.first() )

    emit:
    classifications  = merged_classifications  // channel: [ val(meta), [ kaiju ] , [ kraken ] ]
    kraken           = kraken                  // channel: [ val(meta), [ kraken ] ]
    kaiju            = kaiju                   // channel: [ val(meta), [ kaiju ] ]
    versions         = ch_versions             // channel: [ versions.yml ]
}

