include { BWAMEM2_MEM       } from '../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX     } from '../../modules/nf-core/bwamem2/index/main'
include { BOWTIE2_ALIGN     } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD     } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort/main'

workflow MAP_READS  {

    take:
    reference_reads     // channel: [ val(meta), [ fasta ], [ reads ] ]
    mapper              // val: 'bwamem2' or 'bowtie2'

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    reads       = reference_reads.map{meta, fasta,fastq -> [ meta, fastq ]}
    reference   = reference_reads.map{meta, fasta,fastq -> [ meta, fasta ]}

    if ( mapper == 'bwamem2' ) {
        BWAMEM2_INDEX ( reference )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

        reads_index = reads.join(BWAMEM2_INDEX.out.index, by: [0])
        reads_up    = reads_index.map{meta, reads, index -> [ meta, reads ]}
        index       = reads_index.map{meta, reads, index -> [ meta, index ]}

        BWAMEM2_MEM ( reads_up, index, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
        //no mqc for bwamem2

        ch_bam      = BWAMEM2_MEM.out.bam
    }
    else if ( mapper == 'bowtie2' ) {
        BOWTIE2_BUILD ( reference )
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

        reads_index = reads.join(BOWTIE2_BUILD.out.index, by: [0])
        reads_up    = reads_index.map{meta, reads, index -> [ meta, reads ]}
        index       = reads_index.map{meta, reads, index -> [ meta, index ]}

        BOWTIE2_ALIGN ( reads_up, index, false, false)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bam      = BOWTIE2_ALIGN.out.aligned
        ch_multiqc  = ch_multiqc.mix(BOWTIE2_ALIGN.out.log)

    } else {
        Nextflow.error ("Unknown mapper: ${mapper}")
    }

    SAMTOOLS_SORT ( ch_bam )
    ch_bam      = SAMTOOLS_SORT.out.bam



    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    emit:
    bam      = ch_bam                          // channel: [ val(meta), [ bam ] ]
    ref      = reference                       // channel: [ val(meta), [ fasta ] ]
    mqc      = ch_multiqc                      // channel: [ multiqc ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

