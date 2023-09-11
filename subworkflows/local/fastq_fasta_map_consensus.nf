include { MAP_READS          } from './map_reads'
include { BAM_DEDUPLICATE    } from './bam_deduplicate'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FAIDX     } from '../../modules/nf-core/samtools/faidx/main'
include { BAM_STATS_METRICS  } from './bam_stats_metrics'
include { BAM_CALL_VARIANTS  } from './bam_call_variants'
include { BAM_CALL_CONSENSUS } from './bam_call_consensus'

workflow FASTQ_FASTA_MAP_CONSENSUS {

    take:
    reads                // channel: [ val(meta), [ fastq ] ]
    reference            // channel: [ val(meta), [ fasta ] ]
    mapper               // val: [ bwamem2 | bowtie2 ]
    umi                  // val: [ true | false ]
    deduplicate          // val: [ true | false ]
    variant_caller       // val: [ bcftools | ivar ]
    consensus_caller     // val: [ bcftools | ivar ]
    get_stats            // val: [ true | false ]

    main:

    ch_versions     = Channel.empty()
    ch_multiqc      = Channel.empty()
    ch_dedup_bam    = Channel.empty()

    // We want the meta from the reference channel to be used downstream as this is our varying factor
    // To do this we combine the channels based on sample
    // Extract the reference meta's and reads

    reference
        .map{meta, fasta -> [meta.sample,meta, fasta]}
        .set{ch_reference_tmp}

    reads
        .map{meta, fastq -> [meta.sample,meta, fastq]}
        .set{ch_reads_tmp}

    // Combine the channels based on sample
    ch_reference_tmp
        .combine(ch_reads_tmp, by: [0])
        .map{
            sample, meta_ref, fasta, meta_reads, fastq -> [meta_ref, fasta, fastq]
        }
        .set{ch_reference_reads}

    // mapping of reads using bowtie2 or BWA-MEM2
    MAP_READS ( ch_reference_reads, mapper )

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
    SAMTOOLS_SORT ( ch_dedup_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    ch_dedup_bam_sort = SAMTOOLS_SORT.out.bam

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

    reads_out = ch_reference_reads.map{meta,ref,reads -> [meta,reads] }
    bam_out   = ch_dedup_bam_ref.map{meta,bam,ref -> [meta,bam] }



    emit:
    reads      = reads_out                          // channel: [ val(meta), [ fastq ] ]
    bam        = bam_out                            // channel: [ val(meta), [ bam ] ]
    vcf_filter = ch_vcf_filter                      // channel: [ val(meta), [ vcf ] ]
    vcf        = ch_vcf                             // channel: [ val(meta), [ vcf ] ]
    consensus  = BAM_CALL_CONSENSUS.out.consensus   // channel: [ val(meta), [ fasta ] ]

    mqc      = ch_multiqc                           // channel: [ val(meta), [ csi ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}

