include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { BAM_STATS_SAMTOOLS            } from '../nf-core/bam_stats_samtools/main'

workflow BAM_STATS_METRICS {

    take:
    bam_sort    // channel: [ val(meta), [ bam ] ]
    reference   // channel: [ val(meta), path(fasta) ]
    faidx       // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    SAMTOOLS_INDEX ( bam_sort )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_INDEX.out.bai
        .join(bam_sort)
        .set{ch_bam_sort_bai}

    PICARD_COLLECTMULTIPLEMETRICS ( ch_bam_sort_bai, reference, faidx )
    ch_versions  = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    BAM_STATS_SAMTOOLS ( ch_bam_sort_bai, reference )
    ch_versions  = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.stats.map{it[1]}.ifEmpty{[]})
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.flagstat.map{it[1]}.ifEmpty{[]})
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.idxstats.map{it[1]}.ifEmpty{[]})


    emit:
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    mqc      = ch_multiqc                      // channel: [ multiqc  ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

