include { BAM_VARIANTS_BCFTOOLS } from '../bam_variants_bcftools'
include { BAM_VARIANTS_IVAR     } from '../bam_variants_ivar'
include { VCF_TABIX_STATS       } from '../vcf_tabix_stats'

workflow BAM_CALL_VARIANTS {
    take:
    ch_bam_ref     // channel: [ val(meta), [ bam ], [ fasta ] ]
    variant_caller // value: [ bcftools | ivar ]
    save_stats     // value: [ true | false ]

    main:
    ch_tbi = channel.empty()
    ch_stats = channel.empty()
    ch_multiqc = channel.empty()

    ch_meta_fasta = ch_bam_ref.map { meta, _bam, fasta -> [meta, fasta] }

    if (variant_caller == "bcftools") {
        BAM_VARIANTS_BCFTOOLS(
            ch_bam_ref,
            save_stats,
        )
        ch_vcf = BAM_VARIANTS_BCFTOOLS.out.vcf
        ch_vcf_filter = BAM_VARIANTS_BCFTOOLS.out.vcf_filter
    }
    else if (variant_caller == "ivar") {
        BAM_VARIANTS_IVAR(
            ch_bam_ref,
            save_stats,
        )
        ch_vcf = BAM_VARIANTS_IVAR.out.vcf
        ch_vcf_filter = BAM_VARIANTS_IVAR.out.vcf_filter
        ch_multiqc = ch_multiqc.mix(BAM_VARIANTS_IVAR.out.multiqc)
    }

    if (save_stats) {
        // run stats on all variants not only those that pass the filter
        ch_tabix_in = ch_vcf
            .join(ch_meta_fasta, by: [0])
            .multiMap { meta, vcf, fasta ->
                vcf: [meta, vcf]
                fasta: [meta, fasta]
            }

        VCF_TABIX_STATS(
            ch_tabix_in.vcf,
            [[:], []],
            [[:], []],
            [[:], []],
            [[:], []],
            ch_tabix_in.fasta,
        )
        ch_tbi = VCF_TABIX_STATS.out.index
        ch_stats = VCF_TABIX_STATS.out.stats
        ch_multiqc = ch_multiqc.mix(ch_stats)
    }

    emit:
    vcf        = ch_vcf        // channel: [ val(meta), [ vcf ] ]
    vcf_filter = ch_vcf_filter // channel: [ val(meta), [ vcf ] ]
    tbi        = ch_tbi        // channel: [ val(meta), [ tbi ] ]
    stats      = ch_stats      // channel: [ val(meta), [ stats ] ]
    mqc        = ch_multiqc    // channel: [ val(meta), [ mqc ] ]
}
