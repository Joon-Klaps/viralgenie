include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_1      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_2      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_3      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_4      } from './fastq_fasta_map_consensus.nf'


workflow FASTQ_FASTA_ITERATIVE_CONSENSUS {

    take:
    reference_reads                // channel: [ val(meta), [ reference ], [ fastq ] ]
    repeats                        // val: [ 0 | 1 | 2 | 3 | 4 ]
    intermediate_mapper            // val: [ bwamem2 | bowtie2 ]
    umi                            // val: [ true | false ]
    deduplicate                    // val: [ true | false ]
    call_intermediate_variants     // val: [ true | false ]
    intermediate_variant_caller    // val: [ bcftools | ivar ]
    intermediate_consensus_caller  // val: [ bcftools | ivar ]
    intermediate_mapping_stats     // val: [ true | false ]
    min_mapped_reads               // integer: min_mapped_reads
    min_len                        // integer: min_length
    n_100                          // integer: n_100

    main:
    ch_reference_reads_intermediate = reference_reads
    ch_consensus_allsteps           = Channel.empty()
    ch_multiqc                      = Channel.empty()
    ch_versions                     = Channel.empty()
    ch_bed                          = Channel.empty()
    if (repeats >= 1){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'1', step:"it1", previous_step: meta.step], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_1(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            call_intermediate_variants,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            intermediate_mapping_stats,
            min_mapped_reads,
            min_len,
            n_100
        )

        ch_reference_reads_intermediate = ITERATION_1.out.consensus_reads
        ch_consensus_allsteps           = ch_consensus_allsteps.mix(ITERATION_1.out.consensus)
        ch_multiqc                      = ch_multiqc.mix(ITERATION_1.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_1.out.versions)
        bam                             = ITERATION_1.out.bam
        vcf                             = ITERATION_1.out.vcf
        vcf_filter                      = ITERATION_1.out.vcf_filter
        consensus                       = ITERATION_1.out.consensus
    }
    if (repeats >= 2){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'2', step:"it2", previous_step: meta.step], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_2(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            call_intermediate_variants,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            intermediate_mapping_stats,
            min_mapped_reads,
            min_len,
            n_100
        )

        ch_reference_reads_intermediate = ITERATION_2.out.consensus_reads
        ch_consensus_allsteps           = ch_consensus_allsteps.mix(ITERATION_2.out.consensus)
        ch_multiqc                      = ch_multiqc.mix(ITERATION_2.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_2.out.versions)
        bam                             = ITERATION_2.out.bam
        vcf                             = ITERATION_2.out.vcf
        vcf_filter                      = ITERATION_2.out.vcf_filter
        consensus                       = ITERATION_2.out.consensus
    }
    if (repeats >= 3){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'3', step:"it3", previous_step: meta.step], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_3(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            call_intermediate_variants,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            intermediate_mapping_stats,
            min_mapped_reads,
            min_len,
            n_100
        )

        ch_reference_reads_intermediate = ITERATION_3.out.consensus_reads
        ch_consensus_allsteps           = ch_consensus_allsteps.mix(ITERATION_3.out.consensus)
        ch_multiqc                      = ch_multiqc.mix(ITERATION_3.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_3.out.versions)
        bam                             = ITERATION_3.out.bam
        vcf                             = ITERATION_3.out.vcf
        vcf_filter                      = ITERATION_3.out.vcf_filter
        consensus                       = ITERATION_3.out.consensus
    }
    if (repeats >= 4){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'4', step:"it4", previous_step: meta.step], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_4(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            call_intermediate_variants,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            intermediate_mapping_stats,
            min_mapped_reads,
            min_len,
            n_100
        )

        ch_reference_reads_intermediate = ITERATION_4.out.consensus_reads
        ch_consensus_allsteps           = ch_consensus_allsteps.mix(ITERATION_4.out.consensus)
        ch_multiqc                      = ch_multiqc.mix(ITERATION_4.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_4.out.versions)
        bam                             = ITERATION_4.out.bam
        vcf                             = ITERATION_4.out.vcf
        vcf_filter                      = ITERATION_4.out.vcf_filter
        consensus                       = ITERATION_4.out.consensus
    }

    emit:
    consensus_reads      = ch_reference_reads_intermediate      // channel: [ val(meta), [ fasta ], [ fastq ] ]
    consensus_allsteps   = ch_consensus_allsteps                // channel: [ val(meta), [ fasta ] ]
    bam                  = bam                                  // channel: [ val(meta), [ bam ] ]
    vcf                  = vcf                                  // channel: [ val(meta), [ vcf ] ]
    vcf_filter           = vcf_filter                           // channel: [ val(meta), [ vcf ] ]
    consensus            = consensus                            // channel: [ val(meta), [ fasta ] ]

    mqc                  = ch_multiqc                           // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions                          // channel: [ versions.yml ]
}

