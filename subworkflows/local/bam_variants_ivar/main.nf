// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_ivar.nf

include { IVAR_VARIANTS        } from '../../../modules/nf-core/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF } from '../../../modules/local/ivar_variants_to_vcf'
include { BCFTOOLS_SORT        } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_FILTER      } from '../../../modules/nf-core/bcftools/filter/main'

workflow BAM_VARIANTS_IVAR {
    take:
    ch_bam_fasta // channel: [ val(meta), [ bam ] , [ fasta ]]
    save_stats   // value: [ true | false ]

    main:


    ch_bam     = ch_bam_fasta.map { meta, bam, _fasta -> [meta, bam] }
    ch_fasta   = ch_bam_fasta.map { _meta, _bam, fasta -> [fasta] }
    meta_fasta = ch_bam_fasta.map { meta, _bam, fasta -> [meta, fasta] }

    //
    // Call variants
    //
    IVAR_VARIANTS(
        ch_bam,
        ch_fasta,
        [],
        [],
        save_stats,
    )

    ch_ivar_tsv = IVAR_VARIANTS.out.tsv

    //
    // Convert original iVar output to VCF, zip and index
    //
    ch_ivar_vcf_header = params.ivar_header
        ? file(params.ivar_header, checkIfExists: true)
        : file("${projectDir}/assets/ivar_variants_header_mqc.txt", checkIfExists: true)

    IVAR_VARIANTS_TO_VCF(
        ch_ivar_tsv.join(meta_fasta, by: [0]),
        ch_ivar_vcf_header,
    )


    BCFTOOLS_SORT(
        IVAR_VARIANTS_TO_VCF.out.vcf
    )

    BCFTOOLS_FILTER(
        BCFTOOLS_SORT.out.vcf.map { meta, vcf -> [meta, vcf, []] }
    )

    emit:
    tsv        = ch_ivar_tsv // channel: [ val(meta), [ tsv ] ]
    vcf_orig   = IVAR_VARIANTS_TO_VCF.out.vcf // channel: [ val(meta), [ vcf ] ]
    log_out    = IVAR_VARIANTS_TO_VCF.out.log // channel: [ val(meta), [ log ] ]
    multiqc    = IVAR_VARIANTS_TO_VCF.out.tsv // channel: [ val(meta), [ tsv ] ]
    vcf        = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), [ vcf ] ]
    vcf_filter = BCFTOOLS_FILTER.out.vcf // channel: [ val(meta), [ vcf ] ]
}
