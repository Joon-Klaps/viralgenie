include { MAP_READS                                 } from './map_reads'
include { SAMTOOLS_FAIDX                            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RAW      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP    } from '../../modules/nf-core/samtools/index/main'
include { UMITOOLS_DEDUP                            } from '../../modules/nf-core/umittools/dedup/main'
include { PICARD_DEDUPLICATE                        } from '../../modules/nf-core/umittools/dedup/main'
include { SAMTOOLS_SORT                             } from '../../modules/nf-core/samtools/sort/main'
include { PICARD_COLLECTMULTIPLEMETRICS             } from '../../modules/nf-core/umittools/collectmultiplemetrics/main'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../modules/nf-core/samtools/sort/main'
include { IVAR_CONSENSUS                            } from '../../modules/nf-core/ivar/consensus/main'


workflow FASTQ_FASTA_MAP_CONSENSUS {

    take:
    reads           // channel: [ val(meta), [ fastq ] ]
    reference       // channel: [ val(meta), [ fasta ] ]
    mapper          // val: [ bwamem2 | bowtie2 ]
    umi             // val: [ true | false ]
    deduplicate     // val: [ true | false ]
    get_stats       // val: [ true | false ]

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

    //TODO: Check if all mqc files are in the correct place & modules are called correctly

    MAP_READS ( ch_read_mod, ch_reference_mod, mapper )

    ch_bam       = MAP_READS.out.bam
    ch_versions  =  ch_versions.mix(MAP_READS.out.versions)
    ch_multiqc   =  ch_multiqc.mix(MAP_READS.out.multiqc)

    ch_faidx     = SAMTOOLS_FAIDX ( reference, [[],[]]).faidx

    if (deduplicate) {
        if ( umi ) {
            SAMTOOLS_INDEX_RAW ( MAP_READS.out.bam )
            ch_bam_bai  = SAMTOOLS_INDEX_RAW.out.bai.join(ch_bam, remainder: true)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_RAW.out.versions)

            UMITOOLS_DEDUP ( ch_bam_bai , get_stats)
            ch_dedup_bam  = UMITOOLS_DEDUP.out.bam
            ch_versions   = ch_versions.mix(UMITOOLS_DEDUP.out.versions)


        } else  {
            PICARD_DEDUPLICATE ( MAP_READS.out.bam, reference, ch_faidx )
            ch_dedup_bam      = PICARD_DEDUPLICATE.out.bam
            ch_versions       = ch_versions.mix(PICARD_DEDUPLICATE.out.versions)
        }
    } else {
        ch_dedup_bam = ch_bam
    }

    ch_dedup_sort_bam = SAMTOOLS_SORT( ch_dedup_bam ).sort

    if (get_stats) {
        SAMTOOLS_INDEX_DEDUP ( ch_dedup_sort_bam )
        ch_dedup_sorted_bam_bai  = SAMTOOLS_INDEX_DEDUP.out.bai.join(ch_dedup_sort_bam, remainder: true)
        ch_versions              = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions)

        PICARD_COLLECTMULTIPLEMETRICS ( ch_dedup_sorted_bam_bai, reference, ch_faidx )
        ch_versions  = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
        ch_multiqc   = ch_multiqc.mix(PICARD_COLLECTMULTIPLEMETRICS.out.multiqc)

        BAM_STATS_SAMTOOLS ( ch_dedup_sorted_bam_bai, reference )
        ch_versions  = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.stats)
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.flagstat)
        ch_multiqc   = ch_multiqc.mix(BAM_STATS_SAMTOOLS.out.idxstats)

        //TODO: Check if this is the correct file for multiqc
        ch_multiqc    = ch_multiqc.mix(UMITOOLS_DEDUP.out.log)
    }

    IVAR_CONSENSUS ( ch_dedup_bam, reference.map{it[1]}, get_stats )


    emit:
    //TODO: add if necessary more outputs?
    bam      = ch_dedup_sorted_bam             // channel: [ val(meta), [ bam ] ]
    consensus = IVAR_CONSENSUS.out.consensus  // channel: [ val(meta), [ fasta ] ]

    mqc      = ch_multiqc          // channel: [ val(meta), [ csi ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

