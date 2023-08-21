// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_ivar.nf

include { IVAR_VARIANTS         } from '../../modules/nf-core/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'
include { BCFTOOLS_SORT         } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_FILTER       } from '../../modules/nf-core/bcftools/filter/main'

workflow  {

    take:
    bam         // channel: [ val(meta), [ bam ] ]
    fasta       // channel: [val (meta), [ fasta] ]
    save_stats // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    //
    // Call variants
    //
    IVAR_VARIANTS (
        bam,
        fasta,
        [], // fai never used within ivar
        [], // gff
        save_stats
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())

    // Filter out samples with 0 variants, don't think I want this.
    IVAR_VARIANTS
        .out
        .tsv
     //   .filter { meta, tsv -> WorkflowCommons.getNumLinesInFile(tsv) > 1 }
        .set { ch_ivar_tsv }

    //
    // Convert original iVar output to VCF, zip and index
    //
    IVAR_VARIANTS_TO_VCF (
        ch_ivar_tsv,
        fasta,
        ivar_multiqc_header
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS_TO_VCF.out.versions.first())

    BCFTOOLS_SORT (
        IVAR_VARIANTS_TO_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())


    //
    // TODO: SET modules.config with threshold for variant filtering
    //
    BCFTOOLS_FILTER (
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    emit:
    tsv          = ch_ivar_tsv                     // channel: [ val(meta), [ tsv ] ]

    vcf_orig     = IVAR_VARIANTS_TO_VCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    log_out      = IVAR_VARIANTS_TO_VCF.out.log    // channel: [ val(meta), [ log ] ]
    multiqc      = IVAR_VARIANTS_TO_VCF.out.tsv    // channel: [ val(meta), [ tsv ] ]

    vcf          = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), [ vcf ] ]
    vcf_filter   = BCFTOOLS_FILTER.out.vcf         // channel: [ val(meta), [ vcf ] ]

    versions     = ch_versions                     // channel: [ versions.yml ]
}

// //
// // Function that returns the number of lines in a file
// //
// public static Integer getNumLinesInFile(input_file) {
//     def num_lines = 0
//     input_file.eachLine { line ->
//         num_lines ++
//     }
//     return num_lines
// }

