// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_bcftools.nf
include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM    } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER  } from '../../modules/nf-core/bcftools/filter/main'

workflow BAM_VARIANTS_BCFTOOLS {

    take:
    bam_fasta   // channel: [ val(meta), [ bam ], [ fasta ] ]
    save_stats  // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    ch_bam   = bam_fasta.map{ meta, bam, fasta -> [ meta, bam ] }
    ch_fasta = bam_fasta.map{ meta, bam, fasta -> [ meta, fasta ] }

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP (
        ch_bam.map{ meta, bam_file -> [ meta, bam_file, [] ] },
        ch_fasta,
        save_stats
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // Filter out samples with 0 variants, don't think I wan this?
    BCFTOOLS_MPILEUP
        .out
        .vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .join(ch_fasta)
        .multiMap{
            meta, vcf, tbi, stats, fasta ->
            vcf_tbi : [ meta, vcf, tbi ]
            fasta : [ meta, fasta ]
        }
        .set { bcfnorm_in }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM (
        bcfnorm_in.vcf_tbi,
        bcfnorm_in.fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    //
    // Filter out low quality variants
    //
    BCFTOOLS_FILTER (
        BCFTOOLS_NORM.out.vcf.map{ meta, vcf -> [ meta, vcf, [] ] }
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())


    emit:
    vcf         = BCFTOOLS_NORM.out.vcf    // channel: [ val(meta), [ vcf ] ]
    vcf_filter  = BCFTOOLS_FILTER.out.vcf  // channel: [ val(meta), [ vcf ] ]

    versions    = ch_versions              // channel: [ versions.yml ]
}
