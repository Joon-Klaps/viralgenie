//
// Consensus calling with BCFTools
//

include { HTSLIB_BGZIPTABIX  } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BEDTOOLS_MERGE     } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { MAKE_BED_MASK      } from '../../../modules/local/make_bed_mask/main'

workflow BAM_VCF_CONSENSUS_BCFTOOLS {
    take:
    ch_bam        // channel: [ val(meta), [ bam ] ]
    ch_vcf        // channel: [ val(meta), [ vcf ] ]
    ch_fasta      // channel: [ val(meta), [ fasta ] ]
    mapping_stats // value: [ true | false ]

    main:

    // The VCFs are already bgzipped, so only build the tabix index.
    HTSLIB_BGZIPTABIX(
        ch_vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    ch_bam_vcf_fasta = ch_bam
        .join(ch_vcf, by: [0])
        .join(ch_fasta, by: [0])

    //
    // Create BED file with consensus regions to mask (regions to remove)
    //
    MAKE_BED_MASK(
        ch_bam_vcf_fasta,
        mapping_stats,
    )

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE(
        MAKE_BED_MASK.out.bed
    )

    ch_bed_fasta = BEDTOOLS_MERGE.out.bed
        .join(ch_fasta, by: [0])
        .multiMap {
            meta, bed, fasta ->
                bed : [meta, bed]
                fasta : [fasta]
        }

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA(
        ch_bed_fasta.bed, ch_bed_fasta.fasta
    )

    //
    // Call consensus sequence with BCFTools
    //
    ch_bcftools_in = ch_vcf
        .join(HTSLIB_BGZIPTABIX.out.index, by: [0])
        .join(BEDTOOLS_MASKFASTA.out.fasta, by: [0])
        .map { meta, v, t, f ->
            [meta, v, t, f, []]
        }
    BCFTOOLS_CONSENSUS(
        ch_bcftools_in
    )

    emit:
    consensus = BCFTOOLS_CONSENSUS.out.fasta // channel: [ val(meta), [ fasta ] ]
}
