//
// bin contigs for later reference identification using
//  > Bowtie 2
//  > Metabat2
//

include { BOWTIE2_BUILD                        } from '../../modules/nf-core/bowtie2/build/main'
include { FASTQ_ALIGN_BOWTIE2                  } from '../../subworkflows/nf-core/fastq_align_bowtie2/main'
include { METABAT2_METABAT2                    } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'

workflow FASTA_FASTQ_BOWTIE2_METABAT2 {

    take:
    ch_contigs // channel: [ val(meta), [ fasta ] ]
    ch_reads   // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = channel.empty()

    BOWTIE2_BUILD ( ch_contigs )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    FASTQ_ALIGN_BOWTIE2 ( ch_reads, BOWTIE2_BUILD.out.index, false, true, ch_contigs )
    ch_versions= ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions.first())

    ch_bam_bai = FASTQ_ALIGN_BOWTIE2.out.bam.join(FASTQ_ALIGN_BOWTIE2.out.bai)

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_bam_bai )
    ch_versions= ch_versions.mix( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    ch_contigs_depths= ch_contigs.join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)

    METABAT2_METABAT2 ( ch_contigs )

    emit:
    tooshort        = METABAT2_METABAT2.out.tooshort                    // channel: [ val(meta), [ fasta.gz ] ]
    lowdepth        = METABAT2_METABAT2.out.lowdepth                    // channel: [ val(meta), [ fasta.gz ] ]
    unbinned        = METABAT2_METABAT2.out.unbinned                    // channel: [ val(meta), [ fasta.gz ] ]
    membership      = METABAT2_METABAT2.out.membership                  // channel: [ val(meta), [ tsv.gz ] ]
    bins            = METABAT2_METABAT2.out.fasta                       // channel: [ val(meta), [ fasta.gz ] ]
    contig_depth    = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth    // channel: [ val(meta), [ txt.gz ] ]
    versions        = ch_versions                                       // channel: [ versions.yml ]
}

