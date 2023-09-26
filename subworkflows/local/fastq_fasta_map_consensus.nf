include { MAP_READS                               } from './map_reads'
include { BAM_DEDUPLICATE                         } from './bam_deduplicate'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPPED } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX                          } from '../../modules/nf-core/samtools/faidx/main'
include { BAM_STATS_METRICS                       } from './bam_stats_metrics'
include { BAM_CALL_VARIANTS                       } from './bam_call_variants'
include { BAM_CALL_CONSENSUS                      } from './bam_call_consensus'
include { FASTA_CONTIG_FILTERING                  } from './fasta_contig_filtering'

workflow FASTQ_FASTA_MAP_CONSENSUS {

    take:
    reference_reads      // channel: [ val(meta), [ fasta ], [ fastq ] ]
    mapper               // val: [ bwamem2 | bowtie2 ]
    umi                  // val: [ true | false ]
    deduplicate          // val: [ true | false ]
    variant_caller       // val: [ bcftools | ivar ]
    consensus_caller     // val: [ bcftools | ivar ]
    get_stats            // val: [ true | false ]
    min_len              // integer: min_length
    n_100                // integer: n_100

    main:

    ch_versions     = Channel.empty()
    ch_multiqc      = Channel.empty()
    ch_dedup_bam    = Channel.empty()

    reads_in = reference_reads.map{meta, ref, reads -> [meta,reads] }

    // mapping of reads using bowtie2 or BWA-MEM2
    MAP_READS ( reference_reads, mapper )

    ch_bam       = MAP_READS.out.bam
    ch_reference = MAP_READS.out.ref
    ch_versions  = ch_versions.mix(MAP_READS.out.versions)
    ch_multiqc   = ch_multiqc.mix(MAP_READS.out.mqc)

    SAMTOOLS_FAIDX ( ch_reference, [[],[]])
    ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_bam_fa_fai = ch_bam
        .join(ch_reference, by: [0])
        .join(SAMTOOLS_FAIDX.out.fai, by: [0])

    // deduplicate bam using umitools (if UMI) or picard
    if (deduplicate) {
        BAM_DEDUPLICATE ( ch_bam_fa_fai, umi, get_stats)

        ch_dedup_bam = BAM_DEDUPLICATE.out.bam
        ch_multiqc   = ch_multiqc.mix(BAM_DEDUPLICATE.out.mqc)
        ch_versions  = ch_versions.mix(BAM_DEDUPLICATE.out.versions)

    } else {
        ch_dedup_bam = ch_bam
    }

    // sort bam
    SAMTOOLS_SORT_DEDUPPED ( ch_dedup_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_DEDUPPED.out.versions.first())
    ch_dedup_bam_sort = SAMTOOLS_SORT_DEDUPPED.out.bam

    ch_dedup_bam_ref = ch_dedup_bam_sort.
        join(ch_reference, by: [0])

    // report summary statistics of alignment
    if (get_stats) {
        BAM_STATS_METRICS ( ch_dedup_bam_ref )
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_METRICS.out.mqc)
        ch_versions  = ch_versions.mix(BAM_STATS_METRICS.out.versions)
    }

    // call variants
    ch_vcf        = Channel.empty()
    ch_vcf_filter = Channel.empty()
    ch_tbi        = Channel.empty()

    if (consensus_caller == "bcftools") {
        BAM_CALL_VARIANTS (
            ch_dedup_bam_ref,
            variant_caller,
            get_stats
        )
        ch_versions   = ch_versions.mix(BAM_CALL_VARIANTS.out.versions.first())
        ch_vcf_filter = BAM_CALL_VARIANTS.out.vcf_filter
        ch_vcf        = BAM_CALL_VARIANTS.out.vcf
        ch_tbi        = BAM_CALL_VARIANTS.out.tbi
    }

    // cannot merge ch_vcf_filter as it will not have a meta

    // consensus calling
    BAM_CALL_CONSENSUS (
        ch_dedup_bam_ref,
        ch_vcf_filter,
        consensus_caller,
        get_stats
    )
    ch_versions = ch_versions.mix(BAM_CALL_CONSENSUS.out.versions.first())
    consensus_all   = BAM_CALL_CONSENSUS.out.consensus

    FASTA_CONTIG_FILTERING (
        consensus_all,
        min_len,
        n_100
    )

    consensus_filtered = FASTA_CONTIG_FILTERING.out.contigs

    consensus_reads = consensus_filtered.join(reads_in, by: [0])
    bam_out         = ch_dedup_bam_ref.map{meta,bam,ref -> [meta,bam] }

    emit:
    consensus_reads = consensus_reads                    // channel: [ val(meta), [ fasta ], [ fastq ] ]
    consensus_all   = consensus_all                      // channel: [ val(meta), [ fasta ] ]
    bam             = bam_out                            // channel: [ val(meta), [ bam ] ]
    vcf_filter      = ch_vcf_filter                      // channel: [ val(meta), [ vcf ] ]
    vcf             = ch_vcf                             // channel: [ val(meta), [ vcf ] ]
    consensus       =consensus_filtered                  // channel: [ val(meta), [ fasta ] ]

    mqc             = ch_multiqc                           // channel: [ val(meta), [ csi ] ]
    versions        = ch_versions                          // channel: [ versions.yml ]
}

