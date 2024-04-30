include { KAIJU_KAIJU     as KAIJU_CONTIG      } from '../../modules/nf-core/kaiju/kaiju/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG    } from '../../modules/nf-core/kraken2/kraken2/main'
include { KAIJU_MERGEOUTPUTS                   } from '../../modules/local/kaiju/mergeoutputs/main'
include { EXTRACT_PRECLUSTER                   } from '../../modules/local/extract_precluster/main'

// Classify contigs using kaiju and/or kraken2 and extract their sequences
workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs_reads   // channel: [ val(meta), [ fasta ], [ fastq ] ]
    contig_classifiers // value [kaiju, kraken2]
    ch_kaiju_db        // channel: [ db ]
    ch_kraken2_db      // channel: [ db ]

    main:
    ch_versions = Channel.empty()

    // modify single_end so kaiju & kraken don't crash
    ch_contigs = ch_contigs_reads.map{ meta, fasta,reads -> [meta + [single_end:true, og_single_end:meta.single_end], fasta] }

    kaiju = Channel.empty()
    if ('kaiju' in contig_classifiers){
        KAIJU_CONTIG ( ch_contigs, ch_kaiju_db)
        kaiju       = KAIJU_CONTIG.out.results
        ch_versions = ch_versions.mix( KAIJU_CONTIG.out.versions.first() )
    }

    kraken        = Channel.empty()
    kraken_report = Channel.empty()
    if ('kraken2' in contig_classifiers){
        KRAKEN2_CONTIG ( ch_contigs, ch_kraken2_db, false, true )
        kraken        = KRAKEN2_CONTIG.out.classified_reads_assignment
        kraken_report = KRAKEN2_CONTIG.out.report
        ch_versions   = ch_versions.mix( KRAKEN2_CONTIG.out.versions.first() )
    }

    classifications = Channel.empty()

    if ('kaiju' in contig_classifiers && 'kraken2' in contig_classifiers){
        classifications = kaiju
            .join(kraken, by:[0])
            .join(kraken_report, by:[0])
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kaiju, kraken, kraken_report, contig ->
                kaiju: [meta, kaiju]
                kraken: [meta, kraken, kraken_report]
                contig: [meta, contig]
            }
    } else if ('kaiju' in contig_classifiers){
        classifications = kaiju
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kaiju, contig ->
                kaiju: [meta, kaiju]
                kraken: [meta, [], []]  // empty kraken
                contig: [meta, contig]
            }
    } else if ('kraken2' in contig_classifiers){
        classifications = kraken
            .join(kraken_report, by:[0])
            .join(ch_contigs, by:[0])
            .multiMap{ meta, kraken, kraken_report, contig ->
                kaiju: [meta,[]]    // empty kaiju
                kraken: [meta, kraken, kraken_report]
                contig: [meta, contig]
            }
    }

    EXTRACT_PRECLUSTER ( classifications.kaiju, classifications.kraken, classifications.contig, ch_kaiju_db )
    ch_versions = ch_versions.mix( EXTRACT_PRECLUSTER.out.versions.first() )

    reads = ch_contigs_reads.map{ meta, fasta,reads -> [meta.sample, meta, reads] }
    // modify meta.id to include the taxid & join with reads

    EXTRACT_PRECLUSTER
        .out
        .sequences
        .map { meta, fastas, json ->
            json = WorkflowCommons.getMapFromJson(json)
            return [meta + json, fastas]                                                                // json contains ntaxa
                }
        .transpose()                                                                                    // wide to long
        .map{ meta, fasta ->
            def taxid = fasta.baseName.split("_taxid")[1]                                               // get taxid from fasta file name
            return [meta.sample, meta + [id: "${meta.id}_taxid${taxid}", taxid: "${taxid}"], fasta ]    // [meta.sample, meta, fasta]
        }
        .filter { sample, meta, fasta ->
            params.keep_unclassified || meta.taxid != "U"                                               // filter out unclassified
        }
        .combine(reads, by:[0])                                                                         // reads -> [meta.sample, meta, reads]
        .map{ sample, meta_contig, fasta, meta_reads, reads -> [meta_contig, fasta, reads] }            // select only meta of contigs
        .map{ meta, fasta, reads -> [meta + [single_end:meta.og_single_end], fasta, reads]}             // set original single_end back
        .set{sequences_reads}

    emit:
    sequences_reads  = sequences_reads  // channel: [ [ meta ], [ fasta ], [ fastq ]
    kraken           = kraken           // channel: [ val(meta), [ kraken ] ]
    kaiju            = kaiju            // channel: [ val(meta), [ kaiju ] ]
    versions         = ch_versions      // channel: [ versions.yml ]
}

