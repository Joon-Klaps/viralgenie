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

    bam   = bam_fasta.map{ meta, bam, fasta -> [ meta, bam ] }
    fasta = bam_fasta.map{ meta, bam, fasta -> [ meta, fasta ] }

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP (
        bam.map{ meta, bam_file -> [ meta, bam_file, [] ] },
        fasta,
        save_stats
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // Filter out samples with 0 variants, don't think I wan this?
    BCFTOOLS_MPILEUP
        .out
        .vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .join(fasta)
        // .filter { meta, vcf, tbi, stats -> getNumVariantsFromBCFToolsStats(stats) > 0 }
        .set { ch_vcf_tbi_stats_fasta }

    ch_vcf_tbi_stats_fasta
        .map { meta, vcf, tbi, stats, fasta -> [ meta, vcf, tbi] }
        .set { ch_vcf_tbi }

    // ch_vcf_tbi_stats_fasta
    //     .map { meta, vcf, tbi, stats, fasta -> [ meta, tbi ] }
    //     .set { ch_tbi }

    // ch_vcf_tbi_stats_fasta
    //     .map { meta, vcf, tbi, stats, fasta -> [ meta, stats ] }
    //     .set { ch_stats }

    ch_vcf_tbi_stats_fasta
        .map { meta, vcf, tbi, stats, fasta -> [ meta, fasta ] }
        .set { ch_fasta }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM (
        ch_vcf_tbi,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    //
    // Filter out low quality variants
    //
    BCFTOOLS_FILTER (
        BCFTOOLS_NORM.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())


    emit:
    vcf         = BCFTOOLS_NORM.out.vcf    // channel: [ val(meta), [ vcf ] ]
    vcf_filter  = BCFTOOLS_FILTER.out.vcf  // channel: [ val(meta), [ vcf ] ]

    versions    = ch_versions              // channel: [ versions.yml ]
}

// //
// // Function to get number of variants reported in BCFTools stats file
// //
// public static Integer getNumVariantsFromBCFToolsStats(bcftools_stats) {
//     def num_vars = 0
//     bcftools_stats.eachLine { line ->
//         def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
//         if (matcher) num_vars = matcher[0][1].toInteger()
//     }
//     return num_vars
// }
