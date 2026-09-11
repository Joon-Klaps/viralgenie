//
// Create contigs using
//

include { SPADES                                     } from '../../../modules/nf-core/spades/main'
include { TRINITY                                    } from '../../../modules/nf-core/trinity/main'
include { MEGAHIT                                    } from '../../../modules/nf-core/megahit/main'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_SPADES    } from '../scaffolds_extend_stats'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_TRINITY   } from '../scaffolds_extend_stats'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_MEGAHIT   } from '../scaffolds_extend_stats'
include { FIND_CONCATENATE as CAT_ASSEMBLERS         } from '../../../modules/nf-core/find/concatenate/main'
include { PRINSEQPLUSPLUS as PRINSEQ_CONTIG          } from '../../../modules/nf-core/prinseqplusplus/main'
include { BBMAP_BBNORM                              } from '../../../modules/nf-core/bbmap/bbnorm/main'
include { noContigSamplesToMultiQC                   } from '../utils_nfcore_viralmetagenome_pipeline'


workflow FASTQ_ASSEMBLY {

    take:
    ch_reads        // channel: [ val(meta), [ reads ] ]
    ch_spades_yml   // channel: ['path/to/yml']
    ch_spades_hmm   // channel: ['path/to/hmm']
    normalise_reads // val: [ true | false ] digital normalisation before assembly

    main:
    ch_scaffolds      = channel.empty()
    ch_coverages      = channel.empty()
    ch_multiqc        = channel.empty()
    ch_bad_assemblies = channel.empty()
    assemblers        = params.assemblers ? params.assemblers.split(',').collect{assemblers -> assemblers.trim().toLowerCase() } : []

    // Digital normalisation by k-mer coverage, for the assemblers only.
    ch_reads_assembly = ch_reads
    if (normalise_reads) {
        BBMAP_BBNORM ( ch_reads )
        ch_reads_assembly = BBMAP_BBNORM.out.fastq
    }

    // SPADES
    if ('spades' in assemblers) {
        SPADES(
            ch_reads_assembly.map {meta, reads -> [meta, reads, [], []]},
            ch_spades_yml,
            ch_spades_hmm
            )

        ch_spades_consensus = SPADES.out.scaffolds
            .join(SPADES.out.contigs, remainder:true)
            .map{meta, scaffold, contig -> [meta, scaffold ? scaffold : contig]} // sometimes no scaffold could be created if so take contig

        EXTEND_SPADES( ch_reads, ch_spades_consensus, "spades")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_SPADES.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_SPADES.out.coverages)
        ch_multiqc           = ch_multiqc.mix(EXTEND_SPADES.out.mqc)
    }

    // TRINITY
    if ('trinity' in assemblers) {
        TRINITY(ch_reads_assembly)

        EXTEND_TRINITY( ch_reads, TRINITY.out.transcript_fasta, "trinity")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_TRINITY.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_TRINITY.out.coverages)
        ch_multiqc           = ch_multiqc.mix(EXTEND_TRINITY.out.mqc)
    }

    // MEGAHIT
    if ('megahit' in assemblers) {
        ch_megahit_in = ch_reads_assembly
            .filter {meta, _reads -> meta.single_end }
            .map { meta, reads -> [meta, [reads], []] }
            .mix(
                ch_reads_assembly.filter {meta, _reads -> !meta.single_end }.map { meta, reads -> [meta, [reads[0]], [reads[1]]] }
            )
        MEGAHIT(ch_megahit_in)

        EXTEND_MEGAHIT( ch_reads, MEGAHIT.out.contigs, "megahit")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_MEGAHIT.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_MEGAHIT.out.coverages)
        ch_multiqc           = ch_multiqc.mix(EXTEND_MEGAHIT.out.mqc)
    }

    // ch_scaffolds, go from [[meta,scaffold1],[meta,scaffold2], ...] to [meta,[scaffolds]]
    ch_scaffolds_combined = ch_scaffolds
        .map { meta, scaffold  -> tuple( groupKey(meta, assemblers.size()), scaffold ) }
        .groupTuple(remainder: true)

    ch_coverages_combined = ch_coverages
        .map { meta, coverages  -> tuple( groupKey(meta, assemblers.size()), coverages ) }
        .groupTuple(remainder: true)

    CAT_ASSEMBLERS(ch_scaffolds_combined)
    ch_scaffolds = CAT_ASSEMBLERS.out.file_out

    // Filter out empty scaffolds, might cause certain processes to crash
    ch_scaffolds_branched = ch_scaffolds
        .branch { _meta, scaffolds ->
            pass: scaffolds.countFasta() > 0
            fail: scaffolds.countFasta() == 0
        }

    ch_good_assemblies = ch_scaffolds_branched.pass
    ch_bad_assemblies  = ch_scaffolds_branched.fail

    // Filter low complexity contigs with prinseq++
    if (!params.skip_contig_prinseq){
        ch_prinseq_in = ch_good_assemblies.map{ meta, scaffolds -> [meta, [], scaffolds] }

        PRINSEQ_CONTIG(
            ch_prinseq_in,
        )
        ch_good_assemblies = PRINSEQ_CONTIG.out.good_reads
    }

    ch_good_assemblies_branched = ch_good_assemblies
        .branch { _meta, scaffolds ->
            pass: scaffolds.countFasta() > 0
            fail: scaffolds.countFasta() == 0
        }

    ch_bad_assemblies = ch_bad_assemblies.mix(ch_good_assemblies_branched.fail)
    ch_no_contigs = noContigSamplesToMultiQC(ch_bad_assemblies, assemblers)
        .collectFile(name:'samples_no_contigs_mqc.tsv')
    ch_multiqc = ch_multiqc.mix(ch_no_contigs.ifEmpty([]))

    emit:
    scaffolds            = ch_scaffolds           // channel: [ val(meta), [ scaffolds] ]
    coverages            = ch_coverages_combined  // channel: [ val(meta), [ idxstats* ] ]
    mqc                  = ch_multiqc             // channel: [ val(meta), [ mqc ] ]
    // there are not any MQC files available for spades, trinity and megahit
}
