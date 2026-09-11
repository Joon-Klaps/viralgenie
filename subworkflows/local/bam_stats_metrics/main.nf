include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth/main'
include { BAM_STATS_SAMTOOLS            } from '../../nf-core/bam_stats_samtools/main'
include { CUSTOM_MPILEUP                } from '../../../modules/local/custom_mpileup/main'

workflow BAM_STATS_METRICS {
    take:
    ch_sort_bam_ref // channel: [ val(meta), [ bam ], [ ref ] ]

    main:

    ch_multiqc = channel.empty()

    ch_sort_bam = ch_sort_bam_ref.map { meta, bam, _ref -> [meta, bam] }

    SAMTOOLS_INDEX(ch_sort_bam)

    ch_input_metrics = ch_sort_bam_ref
        .join(SAMTOOLS_INDEX.out.index, by: [0], remainder: true)
        .multiMap { meta, bam, ref, bai ->
            bam_bai: [meta, bam, bai]
            ref: [meta, ref]
            ref_fai: [meta, ref, []]
            bam_bai_bed: [meta, bam, bai, []]
        }

    CUSTOM_MPILEUP(ch_sort_bam_ref)

    PICARD_COLLECTMULTIPLEMETRICS(ch_input_metrics.bam_bai, ch_input_metrics.ref, [[:], []])

    MOSDEPTH(ch_input_metrics.bam_bai_bed, ch_input_metrics.ref, [])
    ch_multiqc  = ch_multiqc.mix(MOSDEPTH.out.global_txt, MOSDEPTH.out.summary_txt)

    BAM_STATS_SAMTOOLS(ch_input_metrics.bam_bai, ch_input_metrics.ref_fai)
    ch_multiqc  = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.stats, BAM_STATS_SAMTOOLS.out.flagstat)

    emit:
    bai      = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ bai ] ]
    mqc      = ch_multiqc // channel: [ multiqc  ]
}
