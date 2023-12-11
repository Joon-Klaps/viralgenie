include { KRAKEN2_BUILD   as KRAKEN2_HOST_BUILD   } from '../../modules/local/kraken2/build/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HOST_REMOVE  } from '../../modules/nf-core/kraken2/kraken2/main'
include { FASTQC          as FASTQC_HOST          } from '../../modules/nf-core/fastqc/main'

workflow FASTQ_KRAKEN_HOST_REMOVE {

    take:
    reads
    kraken2_host_db
    library
    skip_fastqc

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // build kraken2 host db if needed
    if (!kraken2_host_db) {
        KRAKEN2_HOST_BUILD (
            library
        )
        ch_versions     = ch_versions.mix(KRAKEN2_HOST_BUILD.out.versions.first())
        kraken2_host_db = KRAKEN2_HOST_BUILD.out.kraken2_db
    }

    // remove host reads & keep unclassified reads [true, true]
    KRAKEN2_HOST_REMOVE (
        reads,
        kraken2_host_db,
        true,
        true
    )

    ch_versions          = ch_versions.mix(KRAKEN2_HOST_REMOVE.out.versions.first())
    ch_multiqc_files     = ch_multiqc_files.mix( KRAKEN2_HOST_REMOVE.out.report )
    ch_reads_hostremoved = KRAKEN2_HOST_REMOVE.out.unclassified_reads_fastq

    // fastqc
    if (!skip_fastqc) {
        FASTQC_HOST (
            ch_reads_hostremoved
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HOST.out.html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HOST.out.zip)
        ch_versions      = ch_versions.mix(FASTQC_HOST.out.versions.first())
    }

    emit:
    reads_hostremoved = ch_reads_hostremoved            // channel: [ [ meta ], [ fastq ] ]
    mqc               = ch_multiqc_files                // channel: [ multiqc_files ]
    versions          = ch_versions                     // channel: [ versions.yml ]
}

