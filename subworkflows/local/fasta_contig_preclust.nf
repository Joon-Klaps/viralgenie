include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG      } from '../../../modules/nf-core/kraken2/kraken2/main'

workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs         // channel: [ val(meta), [ fasta ] ]
    ch_kraken2_db      // channel: [ val(meta), [ db ] ]

    main:

    ch_versions = Channel.empty()

    KRAKEN2_CONTIG ( reads, ch_kraken2_db, params.kraken2_save_reads, true )
        ch_raw_classifications = ch_raw_classifications.mix(KRAKEN2_KRAKEN2.out.classified_reads_assignment)
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

