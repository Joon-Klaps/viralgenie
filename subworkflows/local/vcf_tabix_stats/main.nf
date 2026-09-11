//
// Run BCFTools tabix and stats commands
// Taken from https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/vcf_tabix_stats.nf
//

include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BCFTOOLS_STATS } from '../../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    ch_vcf     // channel: [ val(meta), [ vcf ] ]
    ch_regions // channel: [ val(meta), [ regions ] ]
    ch_targets // channel: [ val(meta), [ targets ] ]
    ch_samples // channel: [ val(meta), [ samples ] ]
    ch_exons   // channel: [ val(meta), [ exons ] ]
    ch_fasta   // channel: [ val(meta), [ fasta ] ]

    main:

    // The VCFs are already bgzipped, so this only builds the tabix index.
    HTSLIB_BGZIPTABIX(
        ch_vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )
    ch_stats_in = ch_vcf
        .join(HTSLIB_BGZIPTABIX.out.index, by: [0])
        .join(ch_fasta, by: [0])
        .multiMap { meta, vcf, tbi, fasta ->
            vcf_tbi: [meta, vcf, tbi]
            fasta: [meta, fasta]
        }

    BCFTOOLS_STATS(
        ch_stats_in.vcf_tbi,
        ch_regions,
        ch_targets,
        ch_samples,
        ch_exons,
        ch_stats_in.fasta,
    )

    emit:
    index    = HTSLIB_BGZIPTABIX.out.index    // channel: [ val(meta), [ index ] ]
    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]
}
