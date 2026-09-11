include { IVAR_CONSENSUS                                              } from '../../../modules/nf-core/ivar/consensus/main'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CALLED_CONSENSUS } from '../../../modules/local/rename_fasta_header/main'
include { BAM_VCF_CONSENSUS_BCFTOOLS                                  } from '../bam_vcf_consensus_bcftools'

workflow BAM_CALL_CONSENSUS {
    take:
    ch_bam_ref       // channel: [ val(meta), [ bam ], [ fasta ] ]
    ch_vcf           // channel: [ val(meta), [ vcf ] ]
    consensus_caller // value: [ bcftools | ivar ]
    mapping_stats    // value: [ true | false ]

    main:

    ch_bam = ch_bam_ref.map { meta, bam, _fasta -> [meta, bam] }
    ch_fasta = ch_bam_ref.map { meta, _bam, fasta -> [meta, fasta] }

    if (consensus_caller == "bcftools") {
        BAM_VCF_CONSENSUS_BCFTOOLS(
            ch_bam,
            ch_vcf,
            ch_fasta,
            mapping_stats,
        )
        ch_consensus = BAM_VCF_CONSENSUS_BCFTOOLS.out.consensus
    }
    else if (consensus_caller == "ivar") {
        IVAR_CONSENSUS(
            ch_bam,
            ch_fasta.map {_meta, fasta -> fasta },
            mapping_stats,
        )
        ch_consensus = IVAR_CONSENSUS.out.fasta
    }

    RENAME_FASTA_HEADER_CALLED_CONSENSUS(
        ch_consensus,
        consensus_caller,
    )

    emit:
    consensus = RENAME_FASTA_HEADER_CALLED_CONSENSUS.out.fasta // channel: [ val(meta), [ fasta ] ]
}
