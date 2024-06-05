/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

def valid_params = [
    spades_modes     : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
]

def assemblers         = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
def read_classifiers   = params.read_classifiers ? params.read_classifiers.split(',').collect{ it.trim().toLowerCase() } : []
def contig_classifiers = params.precluster_classifiers ? params.precluster_classifiers.split(',').collect{ it.trim().toLowerCase() } : []

def createFileChannel(param) {
    return param ? Channel.fromPath(param, checkIfExists: true).collect() : []
}

def createChannel(dbPath, dbName, skipFlag) {
    return dbPath && skipFlag ? Channel.fromPath(dbPath, checkIfExists: true).map { db -> [[id: dbName], db] } : Channel.empty()
}

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config, params.adapter_fasta, params.contaminants,
    params.spades_yml,params.spades_hmm
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input            ) { ch_input = file(params.input)                                      } else { exit 1, 'Input samplesheet not specified!'                              }

// Optional parameters
ch_adapter_fasta = createFileChannel(params.adapter_fasta)
ch_metadata      = createFileChannel(params.metadata)
ch_contaminants  = createFileChannel(params.contaminants)
ch_spades_yml    = createFileChannel(params.spades_yml)
ch_spades_hmm    = createFileChannel(params.spades_hmm)
ch_constrain_meta = createFileChannel(params.mapping_constrains)

// Databases, we really don't want to stage uncessary databases
ch_ref_pool      = (!params.skip_assembly && !params.skip_polishing) || (!params.skip_consensus_qc && !params.skip_blast_qc)           ? createChannel( params.reference_pool, "reference", true )                                                         : Channel.empty()
ch_kraken2_db    = (!params.skip_assembly && !params.skip_polishing && !params.skip_precluster) || !params.skip_read_classification    ? createChannel( params.kraken2_db, "kraken2", ('kraken2' in read_classifiers || 'kraken2' in contig_classifiers) ) : Channel.empty()
ch_kaiju_db      = (!params.skip_assembly && !params.skip_polishing && !params.skip_precluster) || !params.skip_read_classification    ? createChannel( params.kaiju_db, "kaiju", ('kaiju' in read_classifiers || 'kaiju' in contig_classifiers) )         : Channel.empty()
ch_checkv_db     = !params.skip_consensus_qc                                                                                           ? createChannel( params.checkv_db, "checkv", !params.skip_checkv )                                                  : Channel.empty()
ch_bracken_db    = !params.skip_read_classification                                                                                    ? createChannel( params.bracken_db, "bracken", ('bracken' in read_classifiers) )                                    : Channel.empty()
ch_k2_host       = !params.skip_preprocessing                                                                                          ? createChannel( params.host_k2_db, "k2_host", !params.skip_hostremoval )                                           : Channel.empty()
ch_annotation_db = !params.skip_consensus_qc                                                                                           ? createChannel( params.annotation_db, "annotation", !params.skip_annotation )                                      : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
ch_multiqc_comment_headers            = params.multiqc_comment_headers     ? Channel.fromPath(params.multiqc_comment_headers, checkIfExists:true ) : Channel.empty()
ch_multiqc_custom_table_headers       = params.custom_table_headers        ? Channel.fromPath(params.custom_table_headers, checkIfExists:true ) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL & NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Preprocessing
include { PREPROCESSING_ILLUMINA          } from '../subworkflows/local/preprocessing_illumina'

// metagenomic diversity
include { FASTQ_KRAKEN_KAIJU              } from '../subworkflows/local/fastq_kraken_kaiju'

// Assembly
include { FASTQ_ASSEMBLY    } from '../subworkflows/local/fastq_assembly'
include { noContigSamplesToMultiQC        } from '../modules/local/functions'

// Consensus polishing of genome
include { FASTA_CONTIG_CLUST              } from '../subworkflows/local/fasta_contig_clust'
include { BLAST_MAKEBLASTDB               } from '../modules/nf-core/blast/makeblastdb/main'
include { SEQKIT_REPLACE                  } from '../modules/nf-core/seqkit/replace/main'
include { ALIGN_COLLAPSE_CONTIGS          } from '../subworkflows/local/align_collapse_contigs'
include { UNPACK_DB                       } from '../subworkflows/local/unpack_db'
include { FASTQ_FASTA_ITERATIVE_CONSENSUS } from '../subworkflows/local/fastq_fasta_iterative_consensus'
include { SINGLETON_FILTERING             } from '../subworkflows/local/singleton_filtering'

// Mapping constrains selection
include { FASTQ_FASTA_MASH_SCREEN         } from '../subworkflows/local/fastq_fasta_mash_screen'

// Variant calling
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CONSTRAIN } from '../modules/local/rename_fasta_header'
include { FASTQ_FASTA_MAP_CONSENSUS                            } from '../subworkflows/local/fastq_fasta_map_consensus'

// QC consensus
include { CONSENSUS_QC                    } from '../subworkflows/local/consensus_qc'

// Report generation
include { CUSTOM_MULTIQC_TABLES           } from '../modules/local/custom_multiqc_tables'
include { MULTIQC as MULTIQC_DATAPREP     } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_REPORT       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIRALGENIE {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Importing samplesheet
    ch_reads = Channel.fromSamplesheet(
        'input'
        ).map{
            meta, read1, read2 ->
            single_end = read1 && !read2
            if (single_end) {
                return [meta + [sample: meta.id, single_end: single_end] , [read1]]
            }
            else {
                return [meta + [sample: meta.id, single_end: single_end] , [read1, read2]]
            }
        }

    // Prepare Databases
    ch_db = Channel.empty()
    if ((!params.skip_assembly && !params.skip_polishing) || !params.skip_consensus_qc || !params.skip_read_classification || (!params.skip_preprocessing && !params.skip_hostremoval)){

        ch_db_raw = ch_db.mix(ch_ref_pool,ch_kraken2_db, ch_kaiju_db, ch_checkv_db, ch_bracken_db, ch_k2_host, ch_annotation_db)
        UNPACK_DB (ch_db_raw)

        UNPACK_DB
            .out
            .db
            .branch { meta, unpacked ->
                k2_host: meta.id == 'k2_host'
                    return [ unpacked ]
                reference: meta.id == 'reference'
                    return [ meta, unpacked ]
                checkv: meta.id == 'checkv'
                    return [ unpacked ]
                kraken2: meta.id == 'kraken2'
                    return [ unpacked ]
                bracken: meta.id == 'bracken'
                    return [ unpacked ]
                kaiju: meta.id == 'kaiju'
                    return [ unpacked ]
                annotation: meta.id == 'annotation'
                    return [ meta, unpacked ]
            }
            .set{ch_db}
        ch_versions         = ch_versions.mix(UNPACK_DB.out.versions)

        // transfer to value channels so processes are not just done once
        // '.collect()' is necessary to transform to list so cartesian products are made downstream
        ch_ref_pool_raw     = ch_db.reference.collect{it[1]}.ifEmpty([]).map{it -> [[id: 'reference'], it]}
        ch_kraken2_db       = ch_db.kraken2.collect().ifEmpty([])
        ch_kaiju_db         = ch_db.kaiju.collect().ifEmpty([])
        ch_checkv_db        = ch_db.checkv.collect().ifEmpty([])
        ch_bracken_db       = ch_db.bracken.collect().ifEmpty([])
        ch_k2_host          = ch_db.k2_host.collect().ifEmpty([])
        ch_annotation_db    = ch_db.annotation.collect{it[1]}.ifEmpty([]).map{it -> [[id: 'annotation'], it]}
    }

    // Prepare blast DB
    ch_ref_pool     = Channel.empty()
    ch_blast_refdb  = Channel.empty()
    ch_blast_annodb = Channel.empty()

    if ( (!params.skip_assembly && !params.skip_polishing) || (!params.skip_consensus_qc && !params.skip_blast_qc)){
        ch_blastdb_in = Channel.empty()
        // see issue #56
        SEQKIT_REPLACE (ch_ref_pool_raw)
        ch_versions   = ch_versions.mix(SEQKIT_REPLACE.out.versions)
        ch_ref_pool   = SEQKIT_REPLACE.out.fastx
        ch_blastdb_in = ch_blastdb_in.mix(ch_ref_pool)

        BLAST_MAKEBLASTDB ( ch_blastdb_in )
        BLAST_MAKEBLASTDB
            .out
            .db
            .branch { meta, db ->
                reference: meta.id == 'reference'
                    return [ meta, db ]
                // annotation: meta.id == 'annotation'
                //     return [ meta, db ]
            }.
            set{ch_blastdb_out}
        ch_blast_refdb  = ch_blastdb_out.reference.collect{it[1]}.ifEmpty([]).map{it -> [[id: 'reference'], it]}
        ch_versions     = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    }

    ch_host_trim_reads      = ch_reads
    ch_decomplex_trim_reads = ch_reads
    // preprocessing illumina reads
    if (!params.skip_preprocessing){
        PREPROCESSING_ILLUMINA (
            ch_reads,
            ch_k2_host,
            ch_adapter_fasta,
            ch_contaminants)
        ch_host_trim_reads      = PREPROCESSING_ILLUMINA.out.reads
        ch_decomplex_trim_reads = PREPROCESSING_ILLUMINA.out.reads_decomplexified
        ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.low_reads_mqc.ifEmpty([]))
        ch_versions             = ch_versions.mix(PREPROCESSING_ILLUMINA.out.versions)
    }


    // Determining metagenomic diversity
    if (!params.skip_read_classification) {
        FASTQ_KRAKEN_KAIJU(
            ch_host_trim_reads,
            read_classifiers,
            ch_kraken2_db,
            ch_bracken_db,
            ch_kaiju_db
            )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_KRAKEN_KAIJU.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix(FASTQ_KRAKEN_KAIJU.out.versions)
    }

    // Assembly
    ch_unaligned_raw_contigs   = Channel.empty()
    // Channel for consensus sequences that have been generated across different iteration
    ch_consensus               = Channel.empty()
    // Channel for consensus sequences that have been generated at the LAST iteration
    ch_consensus_results_reads = Channel.empty()
    // Channel for summary table of cluseters to include in mqc report
    ch_clusters_summary        = Channel.empty()

    if (!params.skip_assembly) {
        // run different assemblers and combine contigs
        FASTQ_ASSEMBLY(
            ch_host_trim_reads,
            assemblers,
            ch_spades_yml,
            ch_spades_hmm)

        ch_versions      = ch_versions.mix(FASTQ_ASSEMBLY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ASSEMBLY.out.mqc.collect{it[1]}.ifEmpty([]))

        // Filter out empty scaffolds
        FASTQ_ASSEMBLY
            .out
            .scaffolds
            .branch { meta, scaffolds ->
                pass: scaffolds.countFasta() > 0
                fail: scaffolds.countFasta() == 0
            }
            .set{ch_contigs}

        no_contig_samples = ch_contigs.fail
        noContigSamplesToMultiQC(no_contig_samples, params.assemblers)
            .collectFile(name:'samples_no_contigs_mqc.tsv')
            .set{no_contigs}

        ch_multiqc_files = ch_multiqc_files.mix(no_contigs.ifEmpty([]))


        if (!params.skip_polishing){
            // blast contigs against reference & identify clusters of (contigs & references)
            ch_contigs
                .pass
                .join(ch_host_trim_reads, by: [0], remainder: false)
                .set{ch_contigs_reads}

            FASTA_CONTIG_CLUST (
                ch_contigs_reads,
                ch_blast_refdb,
                ch_ref_pool,
                contig_classifiers,
                ch_kraken2_db,
                ch_kaiju_db
                )
            ch_versions = ch_versions.mix(FASTA_CONTIG_CLUST.out.versions)

            // Split up clusters into singletons and clusters of multiple contigs
            FASTA_CONTIG_CLUST
                .out
                .centroids_members
                .map { meta, centroids, members ->
                    [ meta, centroids, members ]
                }
                .branch { meta, centroids, members ->
                    singletons: meta.cluster_size == 0
                        return [ meta + [step:"singleton"], centroids ]
                    multiple: meta.cluster_size >   0
                        return [ meta + [step:"consensus"], centroids, members ]
                }
                .set{ch_centroids_members}

            ch_clusters_summary    = FASTA_CONTIG_CLUST.out.clusters_summary.collect{it[1]}.ifEmpty([])
            ch_multiqc_files       =  ch_multiqc_files.mix(FASTA_CONTIG_CLUST.out.no_blast_hits_mqc.ifEmpty([]))

            // map clustered contigs & create a single consensus per cluster
            ALIGN_COLLAPSE_CONTIGS (
                ch_centroids_members.multiple
                )
            ch_versions = ch_versions.mix(ALIGN_COLLAPSE_CONTIGS.out.versions)


            SINGLETON_FILTERING (
                ch_centroids_members.singletons,
                params.min_contig_size,
                params.max_n_perc
                )
            ch_versions = ch_versions.mix(SINGLETON_FILTERING.out.versions)

            ALIGN_COLLAPSE_CONTIGS
                .out
                .consensus
                .mix( SINGLETON_FILTERING.out.filtered )
                .set{ ch_consensus }

            ch_unaligned_raw_contigs = ALIGN_COLLAPSE_CONTIGS.out.unaligned_fasta
            ch_unaligned_raw_contigs = ch_unaligned_raw_contigs.mix( SINGLETON_FILTERING.out.filtered )

            // We want the meta from the reference channel to be used downstream as this is our varying factor
            // To do this we combine the channels based on sample
            // Extract the reference meta's and reads
            // Make cartesian product of identified references & reads so all references will be mapped against again.
                ch_decomplex_trim_reads
                    .map{meta, fastq -> [meta.sample,meta, fastq]}
                    .set{ch_reads_tmp}

                ch_consensus
                    .map{meta, fasta -> [meta.sample,meta, fasta]}
                    .set{ch_consensus_tmp}

                ch_consensus_tmp
                    .combine(ch_reads_tmp, by: [0])
                    .map{
                        sample, meta_ref, fasta, meta_reads, fastq -> [meta_ref, fasta, fastq]
                    }
                    .set{ch_consensus_results_reads_intermediate}

            if (!params.skip_iterative_refinement) {
                FASTQ_FASTA_ITERATIVE_CONSENSUS (
                    ch_consensus_results_reads_intermediate,
                    params.iterative_refinement_cycles,
                    params.intermediate_mapper,
                    params.with_umi,
                    params.deduplicate,
                    params.call_intermediate_variants,
                    params.intermediate_variant_caller,
                    params.intermediate_consensus_caller,
                    params.intermediate_mapping_stats,
                    params.min_mapped_reads,
                    params.min_contig_size,
                    params.max_n_perc
                )
                ch_consensus               = ch_consensus.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_allsteps)
                ch_consensus_results_reads = FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_reads
                ch_versions                = ch_versions.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.versions)
                ch_multiqc_files           = ch_multiqc_files.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.mqc.ifEmpty([])) //collect already done in subworkflow
            } else {
                ch_consensus_results_reads = ch_consensus_results_reads_intermediate
            }
        }
    }

    // add last step to it
    ch_consensus_results_reads
        .map{ meta, fasta, fastq ->
            [meta + [step: "variant-calling", iteration:'variant-calling', previous_step: meta.step], fasta, fastq]
            }
        .set{ch_consensus_results_reads}

    if (params.mapping_constrains && !params.skip_variant_calling ) {
        // Importing samplesheet
        Channel.fromSamplesheet('mapping_constrains')
            .map{ meta, sequence ->
                samples = meta.samples == null ? meta.samples: tuple(meta.samples.split(";"))  // Split up samples if meta.samples is not null
                [meta, samples, sequence]
            }
            .transpose(remainder:true)                                                         // Unnest
            .set{ch_mapping_constrains}

        // // Check if the input is a multi-fasta
        // ch_mapping_constrains.map{ meta, samples, sequence -> WorkflowMain.isMultiFasta(sequence, log)}

        //
        ch_decomplex_trim_reads
            .combine( ch_mapping_constrains )
            .filter{ meta_reads, fastq, meta_mapping, mapping_samples, sequence -> mapping_samples == null || mapping_samples == meta_reads.sample}
            .map
                {
                    meta, reads, meta_mapping, samples, sequence_mapping ->
                    id = "${meta.sample}_${meta_mapping.id}-CONSTRAIN"
                    new_meta = meta + meta_mapping + [
                        id: id,
                        cluster_id: "${meta_mapping.id}",
                        step: "constrain",
                        constrain: true,
                        reads: reads,
                        iteration: 'variant-calling',
                        previous_step: 'constrain'
                        ]
                    return [new_meta, sequence_mapping]
                }
            .set{ch_map_seq_anno_combined}

        // Map with both reads and mapping constrains
        ch_map_seq_anno_combined
            .map{ it -> return [it[0], it[1], it[0].reads] }
            .branch{
                meta, fasta, fastq ->
                multiFastaSelection : meta.selection == true
                singleFastaSelection : meta.selection == false
            }
            .set{constrain_consensus_reads}

        // Select the correct reference
        FASTQ_FASTA_MASH_SCREEN (
            constrain_consensus_reads.multiFastaSelection
        )
        ch_versions = ch_versions.mix(FASTQ_FASTA_MASH_SCREEN.out.versions)

        // For QC we keep original sequence to compare to
        ch_unaligned_raw_contigs = ch_unaligned_raw_contigs
            .mix(constrain_consensus_reads.singleFastaSelection.map{meta, fasta, reads -> [meta, fasta]})
            .mix(FASTQ_FASTA_MASH_SCREEN.out.reference_fastq.map{meta, fasta, reads -> [meta, fasta]})

        //Add to the consensus channel, which will be used for variant calling
        ch_consensus_results_reads = ch_consensus_results_reads
            .mix(FASTQ_FASTA_MASH_SCREEN.out.reference_fastq)
            .mix(constrain_consensus_reads.singleFastaSelection)
    }

    // After consensus sequences have been made, we still have to map against it and call variants
    if ( !params.skip_variant_calling ) {

        FASTQ_FASTA_MAP_CONSENSUS(
            ch_consensus_results_reads,
            params.mapper,
            params.with_umi,
            params.deduplicate,
            true,
            params.variant_caller,
            true,
            params.consensus_caller,
            params.mapping_stats,
            params.min_mapped_reads,
            params.min_contig_size,
            params.max_n_perc
        )
        ch_consensus     = ch_consensus.mix(FASTQ_FASTA_MAP_CONSENSUS.out.consensus_all)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTA_MAP_CONSENSUS.out.mqc.ifEmpty([])) // collect already done in subworkflow

    }

    ch_checkv_summary     = Channel.empty()
    ch_quast_summary      = Channel.empty()
    ch_blast_summary      = Channel.empty()
    ch_annotation_summary = Channel.empty()

    if ( !params.skip_consensus_qc  && (!params.skip_assembly || !params.skip_variant_calling) ) {

        CONSENSUS_QC(
            ch_consensus,
            ch_unaligned_raw_contigs,
            ch_checkv_db,
            ch_blast_refdb,
            ch_annotation_db,
            )
        ch_versions           = ch_versions.mix(CONSENSUS_QC.out.versions)
        ch_multiqc_files      = ch_multiqc_files.mix(CONSENSUS_QC.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_checkv_summary     = CONSENSUS_QC.out.checkv_summary.collect{it[1]}.ifEmpty([])
        ch_quast_summary      = CONSENSUS_QC.out.quast_summary.collect{it[1]}.ifEmpty([])
        ch_blast_summary      = CONSENSUS_QC.out.blast_txt.collect{it[1]}.ifEmpty([])
        ch_annotation_summary = CONSENSUS_QC.out.annotation_txt.collect{it[1]}.ifEmpty([])

    }

    MULTIQC_DATAPREP (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
    )

    multiqc_data = MULTIQC_DATAPREP.out.data.ifEmpty([])

    // Prepare MULTIQC custom tables
    CUSTOM_MULTIQC_TABLES (
            ch_clusters_summary.ifEmpty([]),
            ch_metadata,
            ch_checkv_summary.ifEmpty([]),
            ch_quast_summary.ifEmpty([]),
            ch_blast_summary.ifEmpty([]),
            ch_constrain_meta,
            ch_annotation_summary.ifEmpty([]),
            multiqc_data,
            ch_multiqc_comment_headers.ifEmpty([]),
            ch_multiqc_custom_table_headers.ifEmpty([])
            )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQC_TABLES.out.summary_clusters_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQC_TABLES.out.sample_metadata_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQC_TABLES.out.contigs_overview_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQC_TABLES.out.mapping_constrains_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_MULTIQC_TABLES.out.constrains_summary_mqc.ifEmpty([]))
    ch_versions      = ch_versions.mix(CUSTOM_MULTIQC_TABLES.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralgenie.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowViralgenie.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)


    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

    MULTIQC_REPORT (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC_REPORT.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
    if (params.clean_output_on_error) {
        file(params.outdir).deleteDir()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
