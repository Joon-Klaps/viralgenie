include { MAP_READS          } from './map_reads'
include { BAM_DEDUPLICATE    } from './bam_deduplicate'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { BAM_STATS_METRICS  } from './bam_stats_metrics'
include { BAM_CALL_VARIANTS  } from './bam_call_variants'
include { BAM_CALL_CONSENSUS } from './bam_call_consensus'

workflow FASTQ_FASTA_MAP_CONSENSUS {

    take:
    reads            // channel: [ val(meta), [ fastq ] ]
    reference        // channel: [ val(meta), [ fasta ] ]
    mapper           // val: [ bwamem2 | bowtie2 ]
    umi              // val: [ true | false ]
    deduplicate      // val: [ true | false ]
    variant_caller   // val: [ bcftools | ivar ]
    consensus_caller // val: [ bcftools | ivar ]
    get_stats        // val: [ true | false ]

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
        .set{ch_reference_reads}

    // split up in ref & reads
    ch_reference_reads
        .map{
            sample, meta_ref, fasta, meta_reads, fastq -> [meta_ref, fastq]
        }
        .set{ch_read_mod}

    ch_reference_reads
        .map{
            sample, meta_ref, fasta, meta_reads, fastq -> [meta_ref, fasta]
        }
        .set{ch_reference_mod}

    // mapping of reads using bowtie2 or BWA-MEM2
    MAP_READS ( ch_read_mod, ch_reference_mod, mapper )

    ch_bam       = MAP_READS.out.bam
    ch_versions  =  ch_versions.mix(MAP_READS.out.versions)
    ch_multiqc   =  ch_multiqc.mix(MAP_READS.out.mqc)
    ch_faidx     = SAMTOOLS_FAIDX ( reference, [[],[]]).faidx


    // deduplicate bam using umitools (if UMI) or picard
    if (deduplicate) {
        BAM_DEDUPLICATE ( ch_bam, ch_reference_mod, ch_faidx, umi, get_stats)

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

    // report summary statistics of alignment
    if (get_stats) {
        BAM_STATS_METRICS ( ch_dedup_bam_sort, reference, ch_faidx )
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_METRICS.out.mqc)
        ch_versions  = ch_versions.mix(BAM_STATS_METRICS.out.versions)
    }

    // call variants
    ch_vcf = Channel.empty()
    ch_tbi = Channel.empty()

    if (get_stats || consensus_caller == "bcftools") {
        BAM_CALL_VARIANTS( ch_dedup_bam_sort, reference.map{it[1]}, get_stats )
        ch_vcf_filter = BAM_CALL_VARIANTS.out.vcf_filter
        ch_vcf        = BAM_CALL_VARIANTS.out.vcf
        ch_tbi        = BAM_CALL_VARIANTS.out.tbi
    }

    // consensus calling
    BAM_CALL_CONSENSUS (
        ch_dedup_bam_sort,
        ch_vcf_filter,
        ch_reference_mod,
        consensus_caller,
        get_stats
    )


    emit:
    reads      = ch_reference_reads             // channel: [ val(meta), [ fastq ] ]
    bam        = ch_dedup_sorted_bam            // channel: [ val(meta), [ bam ] ]
    vcf_fitler = ch_vcf_filter                  // channel: [ val(meta), [ vcf ] ]
    vcf        = ch_vcf                         // channel: [ val(meta), [ vcf ] ]
    consensus  = BAM_CALL_CONSENSUS.out.consensus   // channel: [ val(meta), [ fasta ] ]

    mqc      = ch_multiqc                       // channel: [ val(meta), [ csi ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}

