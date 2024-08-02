//
// Checks stats on contigs & tries to extend them using SSPACE_BASIC
//
include { QUAST                                } from '../../modules/nf-core/quast/main'
include { SSPACE_BASIC                         } from '../../modules/local/sspace_basic/main'
include { MAP_READS as MAP_READS_CONTIGS       } from './map_reads.nf'
include { SAMTOOLS_INDEX as CONTIG_INDEX       } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS as CONTIG_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'

workflow SCAFFOLDS_EXTEND_STATS {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    scaffolds   // channel: [ val(meta), [ scaffolds ] ]
    name        // value 'spades','trinity','megahit'

    main:
    ch_versions   = Channel.empty()
    ch_scaffolds  = Channel.empty()
    ch_multiqc    = Channel.empty()

    scaffolds
        .filter{meta, contigs -> contigs != null}
        .filter{meta, contigs -> contigs.countFasta() > 0}
        .set{ch_scaffolds}

    // QUAST
    QUAST(ch_scaffolds, [[:],[]], [[:],[]])
    ch_versions = ch_versions.mix(QUAST.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(QUAST.out.tsv.collect{it[1]}.ifEmpty([]))

    // SSPACE_BASIC
    if (!params.skip_sspace_basic){
        ch_scaffolds
            .join(reads)
            .multiMap { meta, scaffolds, reads ->
                reads : [meta, reads]
                scaffolds : [meta, scaffolds]
                settings: [params.read_distance, params.read_distance_sd, params.read_orientation]
                name: name
            }
            .set{ch_sspace_input}

        SSPACE_BASIC(
            ch_sspace_input.reads,
            ch_sspace_input.scaffolds,
            ch_sspace_input.settings,
            ch_sspace_input.name
        )
        ch_versions = ch_versions.mix(SSPACE_BASIC.out.versions.first())

        ch_scaffolds = SSPACE_BASIC.out.scaffolds
    }

    ch_coverages = Channel.empty()
    if (params.perc_reads_contig != 0){
        ch_scaffolds
            .join(reads)
            .set{ch_map_reads_input}

        MAP_READS_CONTIGS(ch_map_reads_input,params.mapper)
        ch_bam      = MAP_READS_CONTIGS.out.bam
        ch_versions = ch_versions.mix(MAP_READS_CONTIGS.out.versions)

        CONTIG_INDEX(ch_bam)
        ch_bam_bai  = ch_bam.join(CONTIG_INDEX.out.bai)
        ch_versions = ch_versions.mix(CONTIG_INDEX.out.versions)

        CONTIG_IDXSTATS(ch_bam_bai)
        ch_coverages = CONTIG_IDXSTATS.out.idxstats
        ch_versions  = ch_versions.mix(CONTIG_IDXSTATS.out.versions)
    }

    emit:
    scaffolds            = ch_scaffolds  // channel: [ val(meta), [ scaffolds] ]
    coverages            = ch_coverages  // channel: [ val(meta), [ idxstats ] ]
    mqc                  = ch_multiqc    // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions   // channel: [ versions.yml ]

}

