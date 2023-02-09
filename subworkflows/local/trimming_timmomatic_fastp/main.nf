//
// A subworkflow that trims the given reads with the given workflow tool
//

include { TRIMMOMATIC                                       } from '../../../modules/nf-core/trimmomatic/main'
include { FASTP                                             } from '../../../modules/nf-core/fastp/main'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_READS   } from '../../../modules/local/multiqc_tsv_from_list'

workflow TRIMMING_TIMMOMATIC_FASTP {

    take:
    reads               // [[meta], [reads]]
    trim_tool           // boolean: 'trimmomatic' or 'fastp'
    adapter_fasta       // fastp: file adapter list

    main:

    trim_reads   = reads
    ch_versions  = Channel.empty()
    trim_log     = Channel.empty()
    trim_summ    = Channel.empty()

    if (trim_tool == 'trimmomatic') {
        TRIMMOMATIC( reads )

        trim_reads  = TRIMMOMATIC.out.trimmed_reads
        trim_log    = TRIMMOMATIC.out.log
        trim_summ   = TRIMMOMATIC.out.summary

        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    } else if  (trim_tool == 'fastp') {
        FASTP( reads, adapter_fasta, false, false )

        trim_log    = FASTP.out.log
        trim_reads  = FASTP.out.reads
        trim_summ   = FASTP.out.json

        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    // //
    // // Filter empty FastQ files after adapter trimming so fastqc will not crash
    // //
    // trim_reads
    //     .join(trim_summ)
    //     .map {
    //         meta, reads, summ ->
    //             pass = PreProcessingTools.getReadsAfterStep(summ, trim_tool) > 0
    //             [ meta, reads, summ, pass ]
    //     }
    //     .set { ch_pass_fail_reads }

    // ch_pass_fail_reads
    //     .map { meta, reads, json, pass -> if (pass) [ meta, reads ] }
    //     .set { trim_reads }

    // ch_pass_fail_reads
    //     .map {
    //         meta, reads, json, pass ->
    //         if (!pass) {
    //             fail_mapped_reads[meta.id] = 0
    //             num_reads = PreProcessingTools.getReadsBeforeStep(summ, trim_tool)
    //             return [ "$meta.id\t$num_reads" ]
    //         }
    //     }
    //     .set { ch_pass_fail_reads }

    // // adding it to the multiqc report
    // MULTIQC_TSV_FAIL_READS (
    //         ch_pass_fail_reads.collect(),
    //         ['Sample', 'Reads before trimming'],
    //         'fail_mapped_reads'
    //     )
    //     .set { fail_reads_multiqc }


    emit:
    reads               = trim_reads            // channel: [ val(meta), [ *.fastq ] ]
    trim_log                                    // channel: [ val(meta), [ *.log ] ]
    //fail_reads_multiqc                          // channel: [ *.tsv ]

    versions            = ch_versions           // channel: [ versions.yml ]
}

