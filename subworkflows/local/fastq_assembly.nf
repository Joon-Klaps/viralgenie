//
// Create contigs using
//

include { SPADES                                     } from '../../modules/nf-core/spades/main'
include { TRINITY                                    } from '../../modules/nf-core/trinity/main'
include { MEGAHIT                                    } from '../../modules/nf-core/megahit/main'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_SPADES    } from './scaffolds_extend_stats.nf'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_TRINITY   } from './scaffolds_extend_stats.nf'
include { SCAFFOLDS_EXTEND_STATS as EXTEND_MEGAHIT   } from './scaffolds_extend_stats.nf'
include { CAT_CAT as CAT_ASSEMBLERS                  } from '../../modules/nf-core/cat/cat/main'
include { PRINSEQPLUSPLUS as PRINSEQ_CONTIG          } from '../../modules/nf-core/prinseqplusplus/main'
include { noContigSamplesToMultiQC                   } from '../../modules/local/functions'


workflow FASTQ_ASSEMBLY {

    take:
    ch_reads        // channel: [ val(meta), [ reads ] ]
    ch_spades_yml   // channel: ['path/to/yml']
    ch_spades_hmm   // channel: ['path/to/hmm']

    main:
    ch_versions    = Channel.empty()
    ch_scaffolds   = Channel.empty()
    ch_coverages   = Channel.empty()
    ch_multiqc     = Channel.empty()
    bad_assemblies = Channel.empty()
    assemblers     = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

    // SPADES
    if ('spades' in assemblers) {
        SPADES(
            ch_reads.map {meta, reads -> [meta, reads, [], []]},
            ch_spades_yml,
            ch_spades_hmm
            )
        ch_versions          = ch_versions.mix(SPADES.out.versions.first())

        SPADES.out.scaffolds
            .join(SPADES.out.contigs, remainder:true)
            .map{meta, scaffold, contig -> [meta, scaffold ? scaffold : contig]} // sometimes no scaffold could be created if so take contig
            .set{spades_consensus}

        EXTEND_SPADES( ch_reads, spades_consensus, "spades")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_SPADES.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_SPADES.out.coverages)
        ch_versions          = ch_versions.mix(EXTEND_SPADES.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_SPADES.out.mqc)
    }

    // TRINITY
    if ('trinity' in assemblers) {
        TRINITY(ch_reads)
        ch_versions          = ch_versions.mix(TRINITY.out.versions.first())

        EXTEND_TRINITY( ch_reads, TRINITY.out.transcript_fasta, "trinity")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_TRINITY.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_TRINITY.out.coverages)
        ch_versions          = ch_versions.mix(EXTEND_TRINITY.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_TRINITY.out.mqc)
    }

    // MEGAHIT
    if ('megahit' in assemblers) {
        MEGAHIT(ch_reads)
        ch_versions          = ch_versions.mix(MEGAHIT.out.versions.first())

        EXTEND_MEGAHIT( ch_reads, MEGAHIT.out.contigs, "megahit")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_MEGAHIT.out.scaffolds)
        ch_coverages         = ch_coverages.mix(EXTEND_MEGAHIT.out.coverages)
        ch_versions          = ch_versions.mix(EXTEND_MEGAHIT.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_MEGAHIT.out.mqc)
    }

    // ch_scaffolds, go from [[meta,scaffold1],[meta,scaffold2], ...] to [meta,[scaffolds]]
    ch_scaffolds
        .map { meta, scaffold  -> tuple( groupKey(meta, assemblers.size()), scaffold ) }
        .groupTuple(remainder: true)
        .set{ch_scaffolds_combined}

    ch_coverages
        .map { meta, coverages  -> tuple( groupKey(meta, assemblers.size()), coverages ) }
        .groupTuple(remainder: true)
        .set{ch_coverages_combined}

    CAT_ASSEMBLERS(ch_scaffolds_combined)
    ch_scaffolds = CAT_ASSEMBLERS.out.file_out
    ch_versions  =  ch_versions.mix(CAT_ASSEMBLERS.out.versions.first())

    // Filter out empty scaffolds, might cause certain processes to crash
    ch_scaffolds
        .branch { meta, scaffolds ->
            pass: scaffolds.countFasta() > 0
            fail: scaffolds.countFasta() == 0
        }
        .set{ch_scaffolds_branched}

    good_assemblies = ch_scaffolds_branched.pass
    bad_assemblies  = ch_scaffolds_branched.fail

    // Filter low complexity contigs with prinseq++
    if (!params.skip_contig_prinseq){
        prinseq_in = good_assemblies.map{ meta, scaffolds -> [meta, [], scaffolds] }

        PRINSEQ_CONTIG(
            prinseq_in,
        )
        ch_versions = ch_versions.mix(PRINSEQ_CONTIG.out.versions.first())
        good_assemblies = PRINSEQ_CONTIG.out.good_reads
    }

    good_assemblies
        .branch { meta, scaffolds ->
            pass: scaffolds.countFasta() > 0
            fail: scaffolds.countFasta() == 0
        }
        .set{good_assemblies_branched}

    bad_assemblies = bad_assemblies.mix(good_assemblies_branched.fail)
    noContigSamplesToMultiQC(bad_assemblies, assemblers)
        .collectFile(name:'samples_no_contigs_mqc.tsv')
        .set{no_contigs}
    ch_multiqc = ch_multiqc.mix(no_contigs.ifEmpty([]))

    emit:
    scaffolds            = ch_scaffolds           // channel: [ val(meta), [ scaffolds] ]
    coverages            = ch_coverages_combined  // channel: [ val(meta), [ idxstats* ] ]
    mqc                  = ch_multiqc             // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions            // channel: [ versions.yml ]
    // there are not any MQC files available for spades, trinity and megahit
}

