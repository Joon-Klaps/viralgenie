//
// Checks stats on contigs & tries to extend them using SSPACE_BASIC
//
include { QUAST                                } from '../../../modules/nf-core/quast/main'
include { SSPACE_BASIC                         } from '../../../modules/local/sspace_basic/main'
include { SAMTOOLS_INDEX as CONTIG_INDEX       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS as CONTIG_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { MAP_READS as MAP_READS_CONTIGS       } from '../map_reads'

workflow SCAFFOLDS_EXTEND_STATS {
    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]
    ch_scaffolds_raw // channel: [ val(meta), [ scaffolds ] ]
    name             // value 'spades','trinity','megahit'

    main:
    ch_scaffolds = channel.empty()
    ch_multiqc = channel.empty()

    ch_scaffolds = ch_scaffolds_raw
        .filter { _meta, contigs -> contigs != null }
        .filter { _meta, contigs -> contigs.countFasta() > 0 }

    // QUAST
    QUAST(ch_scaffolds, [[:], []], [[:], []])
    ch_multiqc = ch_multiqc.mix(QUAST.out.tsv.collect{_meta, tsv -> tsv}.ifEmpty([]))

    // SSPACE_BASIC
    if (!params.skip_sspace_basic) {
        ch_sspace_input = ch_scaffolds
            .join(ch_reads)
            .multiMap { meta, scaffolds, reads ->
                reads: [meta, reads]
                scaffolds: [meta, scaffolds]
                settings: [params.read_distance, params.read_distance_sd, params.read_orientation]
                name: name
            }

        SSPACE_BASIC(
            ch_sspace_input.reads,
            ch_sspace_input.scaffolds,
            ch_sspace_input.settings,
            ch_sspace_input.name,
        )

        ch_scaffolds = SSPACE_BASIC.out.scaffolds
    }

    ch_coverages = channel.empty()
    if (params.perc_reads_contig != 0) {
        ch_map_reads_input = ch_scaffolds.join(ch_reads)

        MAP_READS_CONTIGS(ch_map_reads_input, params.mapper)
        ch_bam = MAP_READS_CONTIGS.out.bam

        CONTIG_INDEX(ch_bam)
        ch_bam_bai = ch_bam.join(CONTIG_INDEX.out.index)

        CONTIG_IDXSTATS(ch_bam_bai)
        ch_coverages = CONTIG_IDXSTATS.out.idxstats
    }

    emit:
    scaffolds = ch_scaffolds // channel: [ val(meta), [ scaffolds] ]
    coverages = ch_coverages // channel: [ val(meta), [ idxstats ] ]
    mqc       = ch_multiqc   // channel: [ val(meta), [ mqc ] ]
}
