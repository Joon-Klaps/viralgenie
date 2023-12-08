include { KAIJU_KAIJU     as KAIJU_CONTIG      } from '../../../modules/nf-core/kaiju/kaiju/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG    } from '../../../modules/nf-core/kraken2/kraken2/main'

workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs         // channel: [ val(meta), [ fasta ] ]
    ch_kaiju_db        // channel: [ val(meta), [ db ] ]
    ch_kraken2_db      // channel: [ val(meta), [ db ] ]

    main:

    ch_versions = Channel.empty()

    KAIJU_CONTIG ( reads, ch_kaiju_db, params.kaiju_save_reads, true )
    ch_classifications_kaiju = KAIJU_CONTIG.out.results
    ch_versions              = ch_versions.mix( KAIJU_CONTIG.out.versions.first() )

    KRAKEN2_CONTIG ( ch_contigs, ch_kraken2_db, false, true )
    ch_classifications_kraken = KRAKEN2_CONTIG.out.classified_reads_assignment
    ch_versions               = ch_versions.mix( KRAKEN2_CONTIG.out.versions.first() )

    ch_classifications = ch_classifications_kaiju.join(ch_classifications_kraken )

    emit:
    classifications_kaiju   = ch_classifications_kaiju      // channel: [ val(meta), [ kaiju ] ]
    classifications_kraken  = ch_classifications_kraken     // channel: [ val(meta), [ kraken ] ]
    classifications         = ch_classifications            // channel: [ val(meta), [ kaiju ] , [ kraken ] ]
    versions                = ch_versions                   // channel: [ versions.yml ]
}

