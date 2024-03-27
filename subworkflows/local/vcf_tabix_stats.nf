//
// Run BCFTools tabix and stats commands
// Taken from https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/vcf_tabix_stats.nf
//

include { TABIX_TABIX    } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    vcf         //    channel: [ val(meta), [ vcf ] ]
    regions     //    channel: [ val(meta), [ regions ] ]
    targets     //    channel: [ val(meta), [ targets ] ]
    samples     //    channel: [ val(meta), [ samples ] ]
    exons       //    channel: [ val(meta), [ exons ] ]
    fasta       //    channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
    vcf_tbi_fasta= vcf
        .join(TABIX_TABIX.out.tbi, by: [0])
        .join(fasta, by: [0])

    vcf_tbi    = vcf_tbi_fasta.map{ meta, vcf, tbi, fasta -> [ meta, vcf, tbi ] }
    meta_fasta = vcf_tbi_fasta.map{ meta, vcf, tbi, fasta -> [ meta, fasta ] }
    BCFTOOLS_STATS (
        vcf_tbi,
        regions,
        targets,
        samples,
        exons,
        meta_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
    tbi      = TABIX_TABIX.out.tbi      // channel: [ val(meta), [ tbi ] ]
    csi      = TABIX_TABIX.out.csi      // channel: [ val(meta), [ csi ] ]

    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    versions = ch_versions              // channel: [ versions.yml ]

}
