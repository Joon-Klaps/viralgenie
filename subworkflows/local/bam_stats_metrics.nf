include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { MOSDEPTH                      } from '../../modules/nf-core/mosdepth/main'
include { BAM_STATS_SAMTOOLS            } from '../nf-core/bam_stats_samtools/main'
include { CUSTOM_MPILEUP                } from '../../modules/local/custom_mpileup/main'

workflow BAM_STATS_METRICS {

    take:
    sort_bam_ref    // channel: [ val(meta), [ bam ], [ ref ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    sort_bam    = sort_bam_ref.map{meta, bam, ref -> [ meta, bam ]}

    SAMTOOLS_INDEX ( sort_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    input_metrics = sort_bam_ref
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .multiMap{
            meta, bam, ref, bai ->
            bam_bai : [meta, bam, bai]
            ref: [meta, ref]
            bam_bai_bed: [meta, bam, bai, []]
        }

    CUSTOM_MPILEUP (sort_bam_ref)
    ch_versions = ch_versions.mix(CUSTOM_MPILEUP.out.versions)

    PICARD_COLLECTMULTIPLEMETRICS ( input_metrics.bam_bai, input_metrics.ref, [[:], []] )
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    MOSDEPTH(input_metrics.bam_bai_bed, input_metrics.ref)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc  = ch_multiqc.mix(MOSDEPTH.out.global_txt)
    ch_multiqc  = ch_multiqc.mix(MOSDEPTH.out.summary_txt)

    BAM_STATS_SAMTOOLS ( input_metrics.bam_bai, input_metrics.ref )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    ch_multiqc  = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.stats)
    ch_multiqc  = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.flagstat)
    // ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.idxstats)

    emit:
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    mqc      = ch_multiqc                      // channel: [ multiqc  ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

