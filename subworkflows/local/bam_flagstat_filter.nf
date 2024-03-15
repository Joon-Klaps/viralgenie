// Take in a bam file and remove those that don't have any reads aligned

include {failedMappedReadsToMultiQC } from '../../modules/local/functions'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT         } from '../../modules/nf-core/samtools/flagstat/main'

workflow BAM_FLAGSTAT_FILTER {

    take:
    ch_bam           // channel: [ val(meta), [ bam ] ]
    min_mapped_reads // integer: min_mapped_reads

    main:

    ch_versions             = Channel.empty()
    ch_fail_mapping_multiqc = Channel.empty()

    SAMTOOLS_INDEX ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai, by: [0])

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_FLAGSTAT
        .out
        .flagstat
        .join(ch_bam, by: [0] )
        .map{ meta, flagstat, bam -> [ meta, bam, WorkflowCommons.getFlagstatMappedReads(flagstat) ] }
        .branch { meta, bam, mapped_reads ->
            pass: mapped_reads > min_mapped_reads
                return [ meta, bam ]
            fail: mapped_reads <= min_mapped_reads
                return [ meta, bam, mapped_reads ]
        }
        .set{ ch_bam_filtered }

    bam_pass = ch_bam_filtered.pass
    bam_fail = ch_bam_filtered.fail

    bam_fail
        .map { meta, bam, mapped_reads ->
            ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.previous_step\t$mapped_reads"]
            }
        .collect()
        .failedMappedReadsToMultiQC()
        .collectFile(name:'failed_mapped_reads_mqc.tsv')
        .set{ch_fail_mapping_multiqc}

    // ch_fail_mapping_multiqc = failedMappedReadsToMultiQC(bam_fail_tsv).collectFile(name:'failed_mapped_reads_mqc.tsv')

    emit:
    bam_pass     = bam_pass                        // channel: [ val(meta), [ bam ] ]
    flagstat     = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    bam_fail_mqc = ch_fail_mapping_multiqc         // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

