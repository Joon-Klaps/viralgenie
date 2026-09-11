include { KAIJU_KAIJU     as KAIJU_CONTIG      } from '../../../modules/nf-core/kaiju/kaiju/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG    } from '../../../modules/nf-core/kraken2/kraken2/main'
include { EXTRACT_PRECLUSTER                   } from '../../../modules/local/extract_precluster/main'
include { getMapFromJson                       } from '../utils_nfcore_viralmetagenome_pipeline'

// Classify contigs using kaiju and/or kraken2 and extract their sequences
workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs_reads   // channel: [ val(meta), [ fasta ], [ fastq ] ]
    contig_classifiers // value:   [ kaiju, kraken2 ]
    ch_kaiju_db        // channel: [ db ]
    ch_kraken2_db      // channel: [ db ]

    main:

    // modify single_end so kaiju & kraken don't crash
    ch_contigs = ch_contigs_reads.map{ meta, fasta, _reads -> [meta + [single_end:true, og_single_end:meta.single_end], fasta] }

    ch_kaiju = channel.empty()
    if ('kaiju' in contig_classifiers){
        KAIJU_CONTIG ( ch_contigs, ch_kaiju_db)
        ch_kaiju    = KAIJU_CONTIG.out.results
    }

    ch_kraken        = channel.empty()
    ch_kraken_report = channel.empty()
    if ('kraken2' in contig_classifiers){
        KRAKEN2_CONTIG ( ch_contigs, ch_kraken2_db, false, true )
        ch_kraken        = KRAKEN2_CONTIG.out.classified_reads_assignment
        ch_kraken_report = KRAKEN2_CONTIG.out.report
    }

    ch_classifications = channel.empty()

    if ('kaiju' in contig_classifiers && 'kraken2' in contig_classifiers){
        ch_classifications = ch_kaiju
            .join(ch_kraken, by:[0])
            .join(ch_kraken_report, by:[0])
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kaiju, kraken, kraken_report, contig ->
                kaiju: [meta, kaiju]
                kraken: [meta, kraken, kraken_report]
                contig: [meta, contig]
            }
    } else if ('kaiju' in contig_classifiers){
        ch_classifications = ch_kaiju
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kaiju, contig ->
                kaiju: [meta, kaiju]
                kraken: [meta, [], []]  // empty kraken
                contig: [meta, contig]
            }
    } else if ('kraken2' in contig_classifiers){
        ch_classifications = ch_kraken
            .join(ch_kraken_report, by:[0])
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kraken, kraken_report, contig ->
                kaiju: [meta,[]]    // empty kaiju
                kraken: [meta, kraken, kraken_report]
                contig: [meta, contig]
            }
    }  else {
        error("No known classifiers found 'kaiju' and 'kraken2' ${contig_classifiers}")
    }

    EXTRACT_PRECLUSTER ( ch_classifications.kaiju, ch_classifications.kraken, ch_classifications.contig, ch_kaiju_db )

    ch_reads = ch_contigs_reads.map{ meta, _fasta, reads -> [meta.sample, meta, reads] }

    // modify meta.id to include the taxid & join with reads
    ch_sequences_reads = EXTRACT_PRECLUSTER
        .out
        .sequences
        .map { meta, fastas, json_file ->
            def json = getMapFromJson(json_file)
            [meta + [ntaxa: json.ntaxa], fastas]                                                        // json contains ONLY ntaxa
        }
        .transpose()                                                                                    // wide to long
        .map{ meta, fasta ->
            def taxid = fasta.baseName.split("_taxid")[1]                                               // get taxid from fasta file name
            return [meta.sample, meta + [id: "${meta.id}_taxid${taxid}", taxid: "${taxid}"], fasta ]    // [meta.sample, meta, fasta]
        }
        .filter { _sample, meta, _fasta ->
            params.keep_unclassified || meta.taxid != "U"                                               // filter out unclassified
        }
        .combine(ch_reads, by:[0])                                                                      // reads -> [meta.sample, meta, reads]
        .map{ _sample, meta_contig, fasta, _meta_reads, reads -> [meta_contig, fasta, reads] }            // select only meta of contigs
        .map{ meta, fasta, reads -> [meta + [single_end:meta.og_single_end], fasta, reads]}             // set original single_end back

    emit:
    contigs_reads  = ch_sequences_reads  // channel: [ [ meta ], [ fasta ], [ fastq ]
    kraken         = ch_kraken           // channel: [ val(meta), [ kraken ] ]
    kaiju          = ch_kaiju            // channel: [ val(meta), [ kaiju ] ]
}
