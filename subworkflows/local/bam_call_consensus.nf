
include { IVAR_CONSENSUS                                                } from '../../modules/nf-core/ivar/consensus/main'
include { BAM_VCF_CONSENSUS_BCFTOOLS                                    } from './bam_vcf_consensus_bcftools.nf'
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CALLED_CONSENSUS   } from '../../modules/local/rename_fasta_header'

workflow BAM_CALL_CONSENSUS {

    take:
    bam_ref          // channel: [ val(meta), [ bam ], [ fasta ] ]
    vcf              // channel: [ val(meta), [ vcf ] ]
    consensus_caller // value: [ bcftools | ivar ]
    mapping_stats        // value: [ true | false ]

    main:

    ch_versions = Channel.empty()
    ch_bam      = bam_ref.map{ meta, bam, fasta -> [ meta, bam ] }
    fasta       = bam_ref.map{ meta, bam, fasta -> [ meta, fasta ] }

    if (consensus_caller == "bcftools"){
        BAM_VCF_CONSENSUS_BCFTOOLS (
            ch_bam,
            vcf,
            fasta,
            mapping_stats
        )
        ch_consensus = BAM_VCF_CONSENSUS_BCFTOOLS.out.consensus
        ch_versions = ch_versions.mix(BAM_VCF_CONSENSUS_BCFTOOLS.out.versions)
    }
    else if (consensus_caller == "ivar"){
        IVAR_CONSENSUS (
            ch_bam,
            fasta.map{it[1]},
            mapping_stats // save mpileup
        )
        ch_consensus = IVAR_CONSENSUS.out.fasta
        ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())
    }

    RENAME_FASTA_HEADER_CALLED_CONSENSUS (
        ch_consensus,
        consensus_caller
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_CALLED_CONSENSUS.out.versions.first())

    emit:
    consensus = RENAME_FASTA_HEADER_CALLED_CONSENSUS.out.fasta   // channel: [ val(meta), [ fasta ] ]
    versions  = ch_versions                                      // channel: [ versions.yml ]
}

