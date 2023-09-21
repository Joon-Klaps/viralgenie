include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_1      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_2      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_3      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_4      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as FINAL_ITERATION  } from './fastq_fasta_map_consensus.nf'

workflow FASTQ_FASTA_ITERATIVE_CONSENSUS {

    take:
    reference_reads                // channel: [ val(meta), [ reference ], [ fastq ] ]
    repeats                        // val: [ 0 | 1 | 2 | 3 | 4 ]
    intermediate_mapper            // val: [ bwamem2 | bowtie2 ]
    final_mapper                   // val: [ bwamem2 | bowtie2 ]
    umi                            // val: [ true | false ]
    deduplicate                    // val: [ true | false ]
    intermediate_variant_caller    // val: [ bcftools | ivar ]
    final_variant_caller           // val: [ bcftools | ivar ]
    intermediate_consensus_caller  // val: [ bcftools | ivar ]
    final_consensus_caller         // val: [ bcftools | ivar ]
    get_intermediate_stats         // val: [ true | false ]
    get_final_stats                // val: [ true | false ]

    main:
    ch_reference_reads_intermediate = reference_reads
    ch_multiqc                      = Channel.empty()
    ch_versions                     = Channel.empty()

    if (repeats > 1){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'1'], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_1(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats
        )

        ch_reference_reads_intermediate = ITERATION_1.out.consensus_reads
        // TODO: Include a check, if coverage is too low, then don't include it in next round, Idea: check if number of ambigous bases is higher than input, if so, fail
        ch_multiqc                      = ch_multiqc.mix(ITERATION_1.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_1.out.versions)
    }
    if (repeats > 2){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'2'], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_2(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats
        )

        ch_reference_reads_intermediate = ITERATION_2.out.consensus_reads
        ch_multiqc                      = ch_multiqc.mix(ITERATION_2.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_2.out.versions)
    }
    if (repeats > 3){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'3'], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        IITERATION_3(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats
        )

        ch_reference_reads_intermediate = ITERATION_3.out.consensus_reads
        ch_multiqc                      = ch_multiqc.mix(ITERATION_3.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_3.out.versions)
    }
    if (repeats > 4){
        ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'4'], fasta, reads]}
            .set{ch_reference_reads_intermediate}

        ITERATION_4(
            ch_reference_reads_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats
        )

        ch_reference_reads_intermediate = ITERATION_4.out.consensus_reads
        ch_multiqc                      = ch_multiqc.mix(ITERATION_4.out.mqc)
        ch_versions                     = ch_versions.mix(ITERATION_4.out.versions)
    }

    ch_reference_reads_intermediate
            .map{meta, fasta, reads -> [meta + [iteration:'final'], fasta, reads]}
            .set{ch_reference_reads_intermediate}

    FINAL_ITERATION(
        ch_reference_reads_intermediate,
        final_mapper,
        umi,
        deduplicate,
        final_variant_caller,
        final_consensus_caller,
        get_final_stats
    )
    ch_multiqc  = ch_multiqc.mix(FINAL_ITERATION.out.mqc)
    ch_versions = ch_versions.mix(FINAL_ITERATION.out.versions)

    emit:
    consensus_reads      = FINAL_ITERATION.out.consensus_reads      // channel: [ val(meta), [ fasta ], [ fastq ] ]
    bam                  = FINAL_ITERATION.out.bam                  // channel: [ val(meta), [ bam ] ]
    vcf                  = FINAL_ITERATION.out.vcf                  // channel: [ val(meta), [ vcf ] ]
    vcf_filter           = FINAL_ITERATION.out.vcf_filter           // channel: [ val(meta), [ vcf ] ]
    consensus            = FINAL_ITERATION.out.consensus            // channel: [ val(meta), [ fasta ] ]

    mqc                  = ch_multiqc                               // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions                              // channel: [ versions.yml ]
}

