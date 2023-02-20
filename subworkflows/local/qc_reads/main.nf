//
// A subworkflow that checks the quality of the given reads with the given workflow tool fastqc or falco
//

include { FASTQC } from '../../../modules/nf-core/fastqc/main'
include { FALCO  } from '../../../modules/nf-core/falco/main'

workflow QC_FASTQC_FALCO {

    take:
    reads //  [ [meta], [ reads ] ]
    qc_tool

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_fastqc_html_report = Channel.empty()

    if (qc_tool == 'fastqc') {
        FASTQC ( reads )
        ch_versions = ch_versions.mix( FASTQC.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip )
        ch_fastqc_html_report = ch_fastqc_html_report.mix( FASTQC.out.html )
    } else if  (qc_tool == 'falco') {
        FALCO ( reads )
        ch_versions = ch_versions.mix( FALCO.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FALCO.out.txt )
        ch_fastqc_html_report = ch_fastqc_html_report.mix( FALCO.out.html )
    }

    emit:
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
    html     = ch_fastqc_html_report
}

