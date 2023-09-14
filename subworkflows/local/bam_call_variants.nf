include { BAM_VARIANTS_BCFTOOLS } from './bam_variants_bcftools.nf'
include { BAM_VARIANTS_IVAR     } from './bam_variants_ivar.nf'
include { VCF_TABIX_STATS       } from './vcf_tabix_stats.nf'

workflow BAM_CALL_VARIANTS {

    take:
    bam_ref         // channel: [ val(meta), [ bam ], [ fasta ] ]
    variant_caller  // value: [ bcftools | ivar ]
    save_stats      // value: [ true | false ]


    main:

    ch_tbi      = Channel.empty()
    ch_csi      = Channel.empty()
    ch_stats    = Channel.empty()
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    bam         = bam_ref.map{ meta, bam, fasta -> [meta, bam ] }
    fasta       = bam_ref.map{ meta, bam, fasta -> fasta }
    meta_fasta  = bam_ref.map{ meta, bam, fasta -> [ meta, fasta ] }

    if (variant_caller == "bcftools"){
        BAM_VARIANTS_BCFTOOLS (
            bam,
            fasta,
            save_stats
        )
        ch_vcf        = BAM_VARIANTS_BCFTOOLS.out.vcf
        ch_vcf_filter = BAM_VARIANTS_BCFTOOLS.out.vcf_filter
        ch_versions   = ch_versions.mix(BAM_VARIANTS_BCFTOOLS.out.versions.first())
    }
    else if (variant_caller == "ivar"){
        BAM_VARIANTS_IVAR (
            bam,
            fasta,
            save_stats
        )
        ch_vcf        = BAM_VARIANTS_IVAR.out.vcf
        ch_vcf_filter = BAM_VARIANTS_IVAR.out.vcf_filter
        ch_versions   = ch_versions.mix(BAM_VARIANTS_IVAR.out.versions.first())
        ch_multiqc    = ch_multiqc.mix(BAM_VARIANTS_IVAR.out.multiqc.map{it[1]}.ifEmpty{[]})
    }

    if (save_stats){
        VCF_TABIX_STATS (
            ch_vcf_filter,
            [[:],[]], // targets
            [[:],[]], // regions
            [[:],[]], // samples
            [[:],[]], // exons
            meta_fasta
        )
        ch_tbi      = VCF_TABIX_STATS.out.tbi
        ch_csi      = VCF_TABIX_STATS.out.csi
        ch_stats    = VCF_TABIX_STATS.out.stats

        ch_versions = ch_versions.mix(VCF_TABIX_STATS.out.versions)
    }

    emit:
    vcf         = ch_vcf         // channel: [ val(meta), [ vcf ] ]
    vcf_filter  = ch_vcf_filter  // channel: [ val(meta), [ vcf ] ]
    tbi         = ch_tbi         // channel: [ val(meta), [ tbi ] ]
    csi         = ch_csi         // channel: [ val(meta), [ csi ] ]
    stats       = ch_stats       // channel: [ val(meta), [ stats ] ]

    mqc         = ch_multiqc     // channel: [ val(meta), [ mqc ] ]
    versions    = ch_versions    // channel: [ versions.yml ]
}

