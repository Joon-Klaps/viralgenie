//
// Read QC, UMI extraction and trimming
//

include { FASTQC as FASTQC_RAW  } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../modules/nf-core/umitools/extract/main'
include { TRIMMOMATIC           } from '../../modules/nf-core/trimmomatic/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getTrimmomaticReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /Input Read Pairs:\s([\d\.]+)/
        def filtered_reads_matcher = line =~ /Dropped Reads:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
}

workflow FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    skip_trimming     // boolean: true/false
    save_trimmed_fail // boolean: true/false
    save_merged       // boolean: true/false
    min_trimmed_reads // integer: > 0

    main:
    ch_versions = Channel.empty()
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC_RAW (
            reads
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip  = FASTQC_RAW.out.zip
        ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    umi_reads = reads
    umi_log   = Channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT (
            reads
        )
        umi_reads   = UMITOOLS_EXTRACT.out.reads
        umi_log     = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1,2]) {
            UMITOOLS_EXTRACT
                .out
                .reads
                .map {
                    meta, reads ->
                        meta.single_end ? [ meta, reads ] : [ meta + [single_end: true], reads[umi_discard_read % 2] ]
                }
                .set { umi_reads }
        }
    }

    trim_reads          = umi_reads
    trim_summ           = Channel.empty()
    trim_unpaired_reads = Channel.empty()
    trim_log            = Channel.empty()
    fastqc_trim_html    = Channel.empty()
    fastqc_trim_zip     = Channel.empty()
    trim_read_count     = Channel.empty()
    if (!skip_trimming) {
        TRIMMOMATIC (
            reads
        )
        trim_summ             = TRIMMOMATIC.out.json
        trim_log              = TRIMMOMATIC.out.log
        trim_unpaired_reads   = TRIMMOMATIC.out.unpaired_reads
        ch_versions           = ch_versions.mix(TRIMMOMATIC.out.versions.first())

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        TRIMMOMATIC
            .out
            .reads
            .join(trim_summ)
            .map { meta, reads, summ -> [ meta, reads, getTrimmomaticReadsAfterFiltering(summ) ] }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { meta, reads, num_reads -> num_reads >= min_trimmed_reads.toInteger() }
            .map { meta, reads, num_reads -> [ meta, reads ] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, reads, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }

        if (!skip_fastqc) {
            FASTQC_TRIM (
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip  = FASTQC_TRIM.out.zip
            ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        }
    }

    emit:
    reads = trim_reads  // channel: [ val(meta), [ reads ] ]

    fastqc_raw_html     // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip      // channel: [ val(meta), [ zip ] ]

    umi_log             // channel: [ val(meta), [ log ] ]

    trim_summ           // channel: [ val(meta), [ summ ] ]
    trim_log            // channel: [ val(meta), [ log ] ]
    trim_unpaired_reads // channel: [ val(meta), [ fastq.gz ] ]
    trim_read_count     // channel: [ val(meta), val(count) ]

    fastqc_trim_html    // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip     // channel: [ val(meta), [ zip ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
