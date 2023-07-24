include { BWAMEM2_MEM    } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX  } from '../../../modules/nf-core/bwamem2/index/main'
include { BOWTIE2_ALIGN  } from '../../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_INDEX  } from '../../../modules/nf-core/bowtie2/index/main'

workflow MAP_READS  {

    take:
    reads     // channel: [ val(meta), [ reads ] ]
    reference // channel: [ val(meta), [ reads ] ]
    mapper       // val: 'bwamem2' or 'bowtie2'

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    if ( mapper == 'bwamem2' ) {
        BWAMEM2_INDEX ( reference )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

        BWAMEM2_MEM ( reads, BWAMEM2_INDEX.out.bwt, true )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

        ch_bam      = BWAMEM2_MEM.out.bam
    }
    else if ( mapper == 'bowtie2' ) {
        BOWTIE2_INDEX ( reference )
        ch_versions = ch_versions.mix(BOWTIE2_INDEX.out.versions.first())

        BOWTIE2_ALIGN ( reads, BOWTIE2_INDEX.out.bt2, true)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bam      = BOWTIE2_ALIGN.out.bam
        ch_multiqc  = ch_multiqc.mix(BOWTIE2_ALIGN.out.log.map{it[1]})

    } else {
        Nextflow.error ("Unknown mapper: ${mapper}")
    }

    emit:
    bam      = ch_bam                          // channel: [ val(meta), [ bam ] ]
    mqc      = ch_multiqc                      // channel: [ val(meta), [ multiqc ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

