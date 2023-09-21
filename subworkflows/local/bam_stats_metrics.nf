include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { BAM_STATS_SAMTOOLS            } from '../nf-core/bam_stats_samtools/main'

workflow BAM_STATS_METRICS {

    take:
    sort_bam_ref    // channel: [ val(meta), [ bam ], [ ref ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    sort_bam    = sort_bam_ref.map{meta, bam, ref -> [ meta, bam ]}
    reference   = sort_bam_ref.map{meta, bam, ref -> [ meta, ref ]}

    SAMTOOLS_INDEX ( sort_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    sort_bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .set{ch_sort_bam_bai}

    PICARD_COLLECTMULTIPLEMETRICS ( ch_sort_bam_bai, reference, [[:], []] )
    ch_versions  = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    BAM_STATS_SAMTOOLS ( ch_sort_bam_bai, reference )
    ch_versions  = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.stats)
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.flagstat)
    ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.idxstats)


    emit:
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    mqc      = ch_multiqc                      // channel: [ multiqc  ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

