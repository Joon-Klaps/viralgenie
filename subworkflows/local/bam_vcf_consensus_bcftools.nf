//
// Consensus calling with BCFTools
//

include { TABIX_TABIX         } from '../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_MERGE      } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA  } from '../../modules/nf-core/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS  } from '../../modules/nf-core/bcftools/consensus/main'
include { MAKE_BED_MASK       } from '../../modules/local/make_bed_mask'
include { RENAME_FASTA_HEADER } from '../../modules/local/rename_fasta_header'

workflow BAM_VCF_CONSENSUS_BCFTOOLS {
    take:
    bam          // channel: [ val(meta), [ bam ] ]
    vcf          // channel: [ val(meta), [ vcf ] ]
    fasta        // channel: [ val(meta), [ fasta ] ]
    get_stats    // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    //
    // Create BED file with consensus regions to mask (regions to remove)
    //
    MAKE_BED_MASK (
        bam.join(vcf, by: [0]),
        fasta,
        get_stats
    )
    ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions.first())

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE (
        MAKE_BED_MASK.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA (
        BEDTOOLS_MERGE.out.bed,
        fasta.map{it[1]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    //
    // Call consensus sequence with BCFTools
    //
    BCFTOOLS_CONSENSUS (
        vcf.join(TABIX_TABIX.out.tbi, by: [0]).join(BEDTOOLS_MASKFASTA.out.fasta, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    //
    // Rename consensus header adding sample name
    //
    RENAME_FASTA_HEADER (
        BCFTOOLS_CONSENSUS.out.fasta
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER.out.versions.first())

    emit:
    consensus        = RENAME_FASTA_HEADER.out.fasta     // channel: [ val(meta), [ fasta ] ]

    versions         = ch_versions                       // channel: [ versions.yml ]
}
