
include { IVAR_CONSENSUS             } from '../../modules/nf-core/ivar/consensus/main'
include { BAM_VCF_CONSENSUS_BCFTOOLS } from './bam_vcf_consensus_bcftools.nf'

workflow BAM_CALL_CONSENSUS {

    take:
    bam              // channel: [ val(meta), [ bam ] ]
    vcf              // channel: [ val(meta), [ vcf ] ]
    fasta            // channel: [val (meta), [ fasta] ]
    consensus_caller // value: [ bcftools | ivar ]
    get_stats        // value: [ true | false ]

    main:

    ch_versions = Channel.empty()

    //TODO: Fix auto completion with correct syntax
    if (consensus_caller == "bcftools"){
        BAM_VCF_CONSENSUS_BCFTOOLS (
            bam,
            vcf,
            fasta,
            get_stats
        )
        ch_consensus = BAM_VCF_CONSENSUS_BCFTOOLS.out.consensus
        ch_versions = ch_versions.mix(BAM_VCF_CONSENSUS_BCFTOOLS.out.versions.first())
    }
    else if (consensus_caller == "ivar"){
        IVAR_CONSENSUS (
            bam,
            fasta,
            get_stats // save mpileup
        )
        ch_consensus = IVAR_CONSENSUS.out.fasta
        ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
    }

    emit:
    consensus = ch_consensus        // channel: [ val(meta), [ fasta ] ]

    versions  = ch_versions          // channel: [ versions.yml ]
}

