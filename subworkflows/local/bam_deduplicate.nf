include { SAMTOOLS_FAIDX     } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP     } from '../../modules/nf-core/umittools/dedup/main'
include { PICARD_DEDUPLICATE } from '../../modules/nf-core/picard/markduplicates/main'

workflow  {

    take:
    bam         // channel: [ val(meta), [ bam ] ]
    reference   // channel: [ val(meta), path(fasta) ]
    faidx       // channel: [ val(meta), path(fasta) ]
    umi         // val: [ true | false ]
    get_stats   // val: [ true | false ]

    main:

    ch_versions = Channel.empty()

    if ( umi ) {
            SAMTOOLS_INDEX( bam )
            ch_bam_bai  = SAMTOOLS_INDEX.out.bai.join(bam, by: [0])
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_RAW.out.versions)

            UMITOOLS_DEDUP ( ch_bam_bai , get_stats)
            ch_dedup_bam  = UMITOOLS_DEDUP.out.bam
            ch_versions   = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
            ch_multiqc    = ch_multiqc.mix(UMITOOLS_DEDUP.out.log)

    } else  {
            PICARD_DEDUPLICATE ( bam, reference, ch_faidx )
            ch_dedup_bam      = PICARD_DEDUPLICATE.out.bam
            ch_versions       = ch_versions.mix(PICARD_DEDUPLICATE.out.versions)
            // no mqc for picard dedup
    }

    emit:
    bam      = ch_dedup_bam           // channel: [ val(meta), [ bam ] ]
    mqc      = ch_multiqc             // channel: [ val(meta), [ multiqc ] ]
    versions = ch_versions            // channel: [ versions.yml ]
}

