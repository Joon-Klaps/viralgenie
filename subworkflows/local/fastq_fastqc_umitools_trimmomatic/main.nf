//
// Read QC, UMI extraction and trimming
//

include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { TRIMMOMATIC           } from '../../../modules/nf-core/trimmomatic/main'


def getTrimmomaticReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /Input Read Pairs:\s([\d\.]+)/
        def filtered_reads_matcher = line =~ /Dropped Reads:\s([\d\.]+)/
        if (total_reads_matcher) {
            total_reads = total_reads_matcher[0][1].toFloat()
        }
        if (filtered_reads_matcher) {
            filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
    }
    return total_reads - filtered_reads
}


workflow FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    skip_trimming     // boolean: true/false
    min_trimmed_reads // integer: > 0

    main:
    ch_fastqc_raw_html = channel.empty()
    ch_fastqc_raw_zip = channel.empty()
    if (!skip_fastqc) {
        FASTQC_RAW(
            ch_reads
        )
        ch_fastqc_raw_html = FASTQC_RAW.out.html
        ch_fastqc_raw_zip = FASTQC_RAW.out.zip
    }

    ch_umi_reads = ch_reads
    ch_umi_log = channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(
            ch_reads
        )
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            ch_umi_reads = UMITOOLS_EXTRACT.out.reads.map { meta, fastq ->
                meta.single_end ? [meta, fastq] : [meta + [single_end: true], fastq[umi_discard_read % 2]]
            }
        }
    }

    ch_trim_reads          = ch_umi_reads
    ch_trim_summ           = channel.empty()
    ch_trim_unpaired_reads = channel.empty()
    ch_trim_log            = channel.empty()
    ch_fastqc_trim_html    = channel.empty()
    ch_fastqc_trim_zip     = channel.empty()
    ch_trim_read_count     = channel.empty()

    if (!skip_trimming) {
        TRIMMOMATIC(
            ch_reads
        )
        ch_trim_summ = TRIMMOMATIC.out.summary
        ch_trim_log = TRIMMOMATIC.out.trim_log
        ch_trim_unpaired_reads = TRIMMOMATIC.out.unpaired_reads

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_num_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
            .join(ch_trim_summ)
            .map { meta, fastq, summ -> [meta, fastq, getTrimmomaticReadsAfterFiltering(summ)] }

        ch_trim_reads = ch_num_trimmed_reads
            .filter { _meta, _fastq, num_reads -> num_reads >= min_trimmed_reads.toInteger() }
            .map { meta, fastq, _num_reads -> [meta, fastq] }

        ch_trim_read_count = ch_num_trimmed_reads.map { meta, _fastq, num_reads -> [meta, num_reads] }

        if (!skip_fastqc) {
            FASTQC_TRIM(
                ch_trim_reads
            )
            ch_fastqc_trim_html = FASTQC_TRIM.out.html
            ch_fastqc_trim_zip = FASTQC_TRIM.out.zip
        }
    }

    emit:
    reads               = ch_trim_reads          // channel: [ val(meta), [ reads ] ]
    fastqc_raw_html     = ch_fastqc_raw_html     // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip      = ch_fastqc_raw_zip      // channel: [ val(meta), [ zip ] ]
    umi_log             = ch_umi_log             // channel: [ val(meta), [ log ] ]
    trim_summ           = ch_trim_summ           // channel: [ val(meta), [ summ ] ]
    trim_log            = ch_trim_log            // channel: [ val(meta), [ log ] ]
    trim_unpaired_reads = ch_trim_unpaired_reads // channel: [ val(meta), [ fastq.gz ] ]
    trim_read_count     = ch_trim_read_count     // channel: [ val(meta), val(count) ]
    fastqc_trim_html    = ch_fastqc_trim_html    // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip     = ch_fastqc_trim_zip     // channel: [ val(meta), [ zip ] ]
}
