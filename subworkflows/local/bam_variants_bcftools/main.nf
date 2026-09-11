// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_bcftools.nf
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM    } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER  } from '../../../modules/nf-core/bcftools/filter/main'

workflow BAM_VARIANTS_BCFTOOLS {
    take:
    ch_bam_fasta // channel: [ val(meta), [ bam ], [ fasta ] ]
    save_stats   // value: [ true | false ]

    main:

    ch_bam = ch_bam_fasta.map { meta, bam, _fasta -> [meta, bam] }
    ch_fasta = ch_bam_fasta.map { meta, _bam, fasta -> [meta, fasta] }
    ch_fasta_fai = ch_bam_fasta.map { meta, _bam, fasta -> [meta, fasta, []] }

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP(
        ch_bam.map { meta, bam_file -> [meta, bam_file, [], []] },
        ch_fasta_fai,
        save_stats,
    )

    // Filter out samples with 0 variants, don't think I wan this?
    ch_bcfnorm_in = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.index)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .join(ch_fasta)
        .multiMap { meta, vcf, tbi, _stats, fasta ->
            vcf_tbi: [meta, vcf, tbi]
            fasta: [meta, fasta]
        }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM(
        ch_bcfnorm_in.vcf_tbi,
        ch_bcfnorm_in.fasta,
    )

    //
    // Filter out low quality variants
    //
    BCFTOOLS_FILTER(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index, by: [0])
    )

    emit:
    vcf        = BCFTOOLS_NORM.out.vcf   // channel: [ val(meta), [ vcf ] ]
    vcf_filter = BCFTOOLS_FILTER.out.vcf // channel: [ val(meta), [ vcf ] ]
}
