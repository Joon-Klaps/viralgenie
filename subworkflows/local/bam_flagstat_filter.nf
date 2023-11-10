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

    ch_versions = Channel.empty()

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
        .branch{ meta, bam, mapped_reads ->
            pass: mapped_reads > min_mapped_reads
                return [ meta + [mapped_reads: mapped_reads], bam ]
            fail: mapped_reads <= min_mapped_reads
                return [ meta + [mapped_reads: mapped_reads], bam ]
        }
        .set{ ch_bam_filtered }


    emit:
    bam_pass = ch_bam_filtered.pass            // channel: [ val(meta), [ bam ] ]
    bam_fail = ch_bam_filtered.fail            // channel: [ val(meta), [ bam ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), [ flagstat ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

