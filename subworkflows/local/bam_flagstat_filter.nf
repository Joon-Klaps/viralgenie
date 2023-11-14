// Take in a bam file and remove those that don't have any reads aligned

include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT  } from '../../modules/nf-core/samtools/flagstat/main'

//
// Function that parses and returns the number of mapped reasds from flagstat files
//
def getFlagstatMappedReads(flagstat_file) {
    def mapped_reads = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped_reads = line.tokenize().first().toInteger()
        }
    }
    return mapped_reads
}

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
        .map{ meta, flagstat, bam -> [ meta, bam, getFlagstatMappedReads(flagstat) ] }
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
            ["$meta.id\t$meta.sample\t$meta.iteration\t$meta.cluster_id\t$mapped_reads"]
            }
        .collect()
        .map {
            tsv_data ->
                def comments = [
                    "id: 'Failed mapped'",
                    "anchor: 'Filtered contigs'",
                    "section_name: 'Minimum mapped reads'",
                    "format: 'tsv'",
                    "description: 'Contigs that did not have more then ${min_mapped_reads} mapped reads were filtered out'",
                    "plot_type: 'table'"
                ]
                def header = ['Id','Sample', 'Iteration','Cluster','Mapped reads']
                return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with other mqc files
        }
        .collectFile(name:'failed_mapped_reads_mqc.tsv')
        .set { ch_fail_mapping_multiqc }

    emit:
    bam_pass     = bam_pass                        // channel: [ val(meta), [ bam ] ]
    flagstat     = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    bam_fail_mqc = ch_fail_mapping_multiqc         // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

