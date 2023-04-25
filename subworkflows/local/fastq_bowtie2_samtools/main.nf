//
// Remove host reads via alignment and export off-target reads
// COPIED & slightly modified from taxprofiler: shortread_hostremoval https://github.com/nf-core/taxprofiler/blob/dev/subworkflows/local/shortread_hostremoval.nf
//

include { BOWTIE2_BUILD             } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN             } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX            } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS            } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_VIEW             } from '../../../modules/nf-core/samtools/view/main'


workflow FASTQ_BOWTIE2_SAMTOOLS {
    take:
    reads       // channel:  [ [ meta ], [ reads ] ]
    reference   //    file:  /path/to/fasta
    index       //    file:  /path/to/index

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( !index ) {
        ch_bowtie2_index = BOWTIE2_BUILD ( [ [], reference ] ).index
        ch_versions      = ch_versions.mix( BOWTIE2_BUILD.out.versions )
    } else {
        ch_bowtie2_index = index.first()
    }

    // Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
    BOWTIE2_ALIGN ( reads, ch_bowtie2_index, true, true)
    ch_versions      = ch_versions.mix( BOWTIE2_ALIGN.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( BOWTIE2_ALIGN.out.log )

    // Indexing whole BAM for host removal statistics
    SAMTOOLS_INDEX ( BOWTIE2_ALIGN.out.bam )
    ch_versions      = ch_versions.mix( SAMTOOLS_INDEX.out.versions.first() )

    bam_bai = BOWTIE2_ALIGN.out.bam
        .join(SAMTOOLS_INDEX.out.bai, remainder: true)

    SAMTOOLS_STATS ( bam_bai, reference )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_STATS.out.stats )

    emit:
    stats    = SAMTOOLS_STATS.out.stats
    reads    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions               // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
