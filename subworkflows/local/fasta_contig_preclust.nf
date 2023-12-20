include { KAIJU_KAIJU     as KAIJU_CONTIG      } from '../../modules/nf-core/kaiju/kaiju/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIG    } from '../../modules/nf-core/kraken2/kraken2/main'
include { KAIJU_MERGEOUTPUTS                   } from '../../modules/local/kaiju/mergeoutputs/main'
include { EXTRACT_PRECLUSTER                   } from '../../modules/local/extract_precluster/main'

// Classify contigs using kaiju and/or kraken2 and extract their sequences
workflow FASTA_CONTIG_PRECLUST {

    take:
    ch_contigs_reads   // channel: [ val(meta), [ fasta ], [ fastq ] ]
    ch_kaiju_db        // channel: [ val(meta), [ db ] ]
    ch_kraken2_db      // channel: [ val(meta), [ db ] ]

    main:
    ch_versions = Channel.empty()

    // modify single_end so kaiju & kraken don't crash
    ch_contigs = ch_contigs_reads.map{ meta, fasta,reads -> [meta + [single_end:true, og_single_end:meta.single_end], fasta] }

    kaiju = Channel.empty()
    if (!params.skip_kaiju){
        KAIJU_CONTIG ( ch_contigs, ch_kaiju_db)
        kaiju       = KAIJU_CONTIG.out.results
        ch_versions = ch_versions.mix( KAIJU_CONTIG.out.versions.first() )
    }

    kraken = Channel.empty()
    if (!params.skip_kraken2){
        KRAKEN2_CONTIG ( ch_contigs, ch_kraken2_db, false, true )
        kraken       = KRAKEN2_CONTIG.out.classified_reads_assignment
        ch_versions  = ch_versions.mix( KRAKEN2_CONTIG.out.versions.first() )
    }

    classifications = kaiju.join(kraken, remainder:true)

    if (!(params.skip_kaiju || params.skip_kraken2)){
        KAIJU_MERGEOUTPUTS ( classifications, ch_kaiju_db )
        classifications = KAIJU_MERGEOUTPUTS.out.merged
        ch_versions     = ch_versions.mix( KAIJU_MERGEOUTPUTS.out.versions.first() )
    } else if (!params.skip_kaiju){
        classifications = kaiju
    } else if (!params.skip_kraken2){
        classifications = kraken
    }

    clas_contig = classifications.join(ch_contigs, by:[0])

    EXTRACT_PRECLUSTER ( clas_contig )
    ch_versions = ch_versions.mix( EXTRACT_PRECLUSTER.out.versions.first() )

    reads = ch_contigs_reads.map{ meta, fasta,reads -> [meta.sample, meta, reads] }
    // modify meta.id to include the taxid & join with reads

    EXTRACT_PRECLUSTER
        .out
        .sequences
        .map { meta, fastas, json ->
            json = WorkflowCommons.getMapFromJson(json)
            return [meta + json, fastas]                                                        // json contains ntaxa
                }
        .transpose()                                                                            // wide to long
        .map{ meta, fasta ->
            def taxid = fasta.baseName.split("_taxid")[1]                                       // get taxid from fasta file name
            return [meta.sample, meta + [id: "${meta.id}_taxid${taxid}"], fasta ]               // [meta.sample, meta, fasta]
        }
        .combine(reads, by:[0])                                                                 // reads -> [meta.sample, meta, reads]
        .map{ sample, meta_contig, fasta, meta_reads, reads -> [meta_contig, fasta, reads] }    // select only meta of contigs
        .map{ meta, fasta, reads -> [meta + [single_end:meta.og_single_end], fasta, reads]}     // set original single_end back
        .set{sequences_reads}

    classifications = classifications.map{ meta, txt -> [meta + [single_end:meta.og_single_end], txt] }

    emit:
    sequences_reads  = sequences_reads  // channel: [ [ meta ], [ fasta ], [ fastq ]
    classifications  = classifications  // channel: [ val(meta), [ kaiju ] , [ kraken ] ]
    kraken           = kraken           // channel: [ val(meta), [ kraken ] ]
    kaiju            = kaiju            // channel: [ val(meta), [ kaiju ] ]
    versions         = ch_versions      // channel: [ versions.yml ]
}

