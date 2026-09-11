include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP        } from '../../../modules/nf-core/umitools/dedup/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'

workflow BAM_DEDUPLICATE {
    take:
    ch_bam_ref_fai // channel: [ val(meta), [ bam ], [ fasta ], [ fai ] ]
    umi            // val: [ true | false ]
    mapping_stats  // val: [ true | false ]

    main:

    ch_multiqc = channel.empty()

    ch_bam = ch_bam_ref_fai.map { meta, bam, _fasta, _fai -> [meta, bam] }
    ch_reference_fai = ch_bam_ref_fai.map { meta, _bam, fasta, fai -> [meta, fasta, fai] }

    if (umi && ['mapping', 'both'].contains(params.umi_deduplicate)) {
        SAMTOOLS_INDEX(ch_bam)
        ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.index, by: [0])

        UMITOOLS_DEDUP(ch_bam_bai, mapping_stats)
        ch_dedup_bam = UMITOOLS_DEDUP.out.bam
        if (mapping_stats) {
            ch_multiqc = ch_multiqc.mix(UMITOOLS_DEDUP.out.log)
        }
    }
    else {
        PICARD_MARKDUPLICATES(ch_bam, ch_reference_fai)
        ch_dedup_bam = PICARD_MARKDUPLICATES.out.bam
        if (mapping_stats) {
            ch_multiqc = ch_multiqc.mix(PICARD_MARKDUPLICATES.out.metrics)
        }
    }

    emit:
    bam      = ch_dedup_bam // channel: [ val(meta), [ bam ] ]
    mqc      = ch_multiqc   // channel: [ multiqc ]
}
