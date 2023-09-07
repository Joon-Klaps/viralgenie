include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_1      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_2      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_3      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_4      } from './fastq_fasta_map_consensus.nf'
include { FASTQ_FASTA_MAP_CONSENSUS as ITERATION_FINAL  } from './fastq_fasta_map_consensus.nf'

workflow FASTQ_FASTA_ITERATIVE_CONSENSUS {

    take:
    reads                          // channel: [ val(meta), [ fastq ] ]
    reference                      // channel: [ val(meta), [ fasta ] ]
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
    ivar_header                    // path: [ header ]

    main:

    ch_reference_intermediate = reference
    ch_multiqc                = Channel.empty()
    ch_versions               = Channel.empty()

    if (repeats > 1){
        ch_reference_intermediate
            .map{meta, fasta -> [meta + [iteration:'1'], fasta]}
            .set{ch_reference_intermediate}

        ch_reference_intermediate.view()
        ITERATION_1(
            reads,
            ch_reference_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats,
            ivar_header
        )

        ch_reference_intermediate = ITERATION_1.out.consensus
        ch_multiqc                = ch_multiqc.mix(ITERATION_1.out.mqc)
        ch_versions               = ch_versions.mix(ITERATION_1.out.versions)
    }
    if (repeats > 2){
        ch_reference_intermediate
            .map{meta, fasta -> [meta + [iteration:'2'], fasta]}
            .set{ch_reference_intermediate}

        ITERATION_2(
            reads,
            ch_reference_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats,
            ivar_header
        )

        ch_reference_intermediate = ITERATION_2.out.consensus
        ch_multiqc                = ch_multiqc.mix(ITERATION_2.out.mqc)
    }
    if (repeats > 3){
        ch_reference_intermediate
            .map{meta, fasta -> [meta + [iteration:'3'], fasta]}
            .set{ch_reference_intermediate}

        IITERATION_3(
            reads,
            ch_reference_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats,
            ivar_header
        )

        ch_reference_intermediate = ITERATION_3.out.consensus
        ch_multiqc                = ch_multiqc.mix(ITERATION_3.out.mqc)
    }
    if (repeats > 4){
        ch_reference_intermediate
            .map{meta, fasta -> [meta + [iteration:'4'], fasta]}
            .set{ch_reference_intermediate}

        ITERATION_4(
            reads,
            ch_reference_intermediate,
            intermediate_mapper,
            umi,
            deduplicate,
            intermediate_variant_caller,
            intermediate_consensus_caller,
            get_intermediate_stats,
            ivar_header
        )

        ch_reference_intermediate = ITERATION_4.out.consensus
        ch_multiqc                = ch_multiqc.mix(ITERATION_4.out.mqc)
    }

    ch_reference_intermediate
            .map{meta, fasta -> [meta + [iteration:'final'], fasta]}
            .set{ch_reference_intermediate}

    ITERATION_FINAL(
        reads,
        ch_reference_intermediate,
        final_mapper,
        umi,
        deduplicate,
        final_variant_caller,
        final_consensus_caller,
        get_final_stats,
        ivar_header
    )

    ch_multiqc = ch_multiqc.mix(ITERATION_FINAL.out.mqc)
    ch_versions = ch_versions.mix(ITERATION_FINAL.out.versions)


    emit:
    reads      = ITERATION_FINAL.out.reads      // channel: [ val(meta), [ fastq ] ]
    bam        = ITERATION_FINAL.out.bam        // channel: [ val(meta), [ bam ] ]
    vcf        = ITERATION_FINAL.out.vcf        // channel: [ val(meta), [ vcf ] ]
    vcf_filter = ITERATION_FINAL.out.vcf_filter // channel: [ val(meta), [ vcf ] ]
    consensus  = ITERATION_FINAL.out.consensus  // channel: [ val(meta), [ fasta ] ]

    mqc        = ch_multiqc                     // channel: [ val(meta), [ mqc ] ]
    versions   = ch_versions                    // channel: [ versions.yml ]
}

