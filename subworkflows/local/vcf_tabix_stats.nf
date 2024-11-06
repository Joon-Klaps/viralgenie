//
// Run BCFTools tabix and stats commands
// Taken from https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/vcf_tabix_stats.nf
//

include { TABIX_TABIX    } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    ch_vcf      // channel: [ val(meta), [ vcf ] ]
    ch_regions  // channel: [ val(meta), [ regions ] ]
    ch_targets  // channel: [ val(meta), [ targets ] ]
    ch_samples  // channel: [ val(meta), [ samples ] ]
    ch_exons    // channel: [ val(meta), [ exons ] ]
    ch_fasta    // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        ch_vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
    ch_vcf
        .join(TABIX_TABIX.out.tbi, by: [0])
        .join(ch_fasta, by: [0])
        .multiMap{meta, vcf, tbi, fasta ->
            vcf_tbi : [ meta, vcf, tbi]
            fasta : [ meta, fasta ]
        }
        .set{stats_in}

    BCFTOOLS_STATS (
        stats_in.vcf_tbi,
        ch_regions,
        ch_targets,
        ch_samples,
        ch_exons,
        stats_in.meta_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
    tbi      = TABIX_TABIX.out.tbi      // channel: [ val(meta), [ tbi ] ]
    csi      = TABIX_TABIX.out.csi      // channel: [ val(meta), [ csi ] ]

    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    versions = ch_versions              // channel: [ versions.yml ]

}
