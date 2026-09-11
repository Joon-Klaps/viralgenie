/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VIRALMETAGENOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL & NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// functions
include { samplesheetToList               } from 'plugin/nf-schema'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
include { createFileChannel               } from '../subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
include { createChannel                   } from '../subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
include { noContigSamplesToMultiQC        } from '../subworkflows/local/utils_nfcore_viralmetagenome_pipeline'
include { getLengthAndAmbigous            } from '../subworkflows/local/utils_nfcore_viralmetagenome_pipeline'

// Preprocessing
include { PREPROCESSING_ILLUMINA          } from '../subworkflows/local/preprocessing_illumina'
// metagenomic diversity
include { FASTQ_KRAKEN_KAIJU              } from '../subworkflows/local/fastq_kraken_kaiju'
// Assembly
include { FASTQ_ASSEMBLY                  } from '../subworkflows/local/fastq_assembly'
// Consensus polishing of genome
include { FASTA_CONTIG_CLUST              } from '../subworkflows/local/fasta_contig_clust'
include { BLAST_MAKEBLASTDB               } from '../modules/nf-core/blast/makeblastdb/main'
include { ALIGN_COLLAPSE_CONTIGS          } from '../subworkflows/local/align_collapse_contigs'
include { UNPACK_DB                       } from '../subworkflows/local/unpack_db'
include { FASTQ_FASTA_ITERATIVE_CONSENSUS } from '../subworkflows/local/fastq_fasta_iterative_consensus'
include { SINGLETON_FILTERING             } from '../subworkflows/local/singleton_filtering'
// Mapping constraints selection
include { FASTQ_FASTA_MASH_SCREEN         } from '../subworkflows/local/fastq_fasta_mash_screen'
// QC consensus
include { CONSENSUS_QC                    } from '../subworkflows/local/consensus_qc'
// Report generation
include { CUSTOM_MULTIQC                  } from '../modules/local/custom_multiqc'
// Variant calling
include { FASTQ_FASTA_MAP_CONSENSUS       } from '../subworkflows/local/fastq_fasta_map_consensus'
include { VCF_ANNOTATE                    } from '../subworkflows/local/vcf_annotate'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRALMETAGENOME {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir
    keep_unmapped // boolean: keep unmapped reads in downstream alignments
    normalise_reads // boolean: digitally normalise reads before assembly
    use_host_filtered_reads // boolean: prefer host-filtered reads for downstream mapping & polishing steps

    main:

    ch_multiqc_files = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PARAMETER INITIALIZATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def read_classifiers   = params.read_classifiers ? params.read_classifiers.split(',').collect{classifiers -> classifiers.trim().toLowerCase() } : []
    def contig_classifiers = params.precluster_classifiers ? params.precluster_classifiers.split(',').collect{classifiers -> classifiers.trim().toLowerCase() } : []
    // Optional parameters
    ch_blacklist       = createFileChannel(params.blacklist)
    ch_metadata        = createFileChannel(params.metadata)
    ch_contaminants    = createFileChannel(params.contaminants)
    ch_spades_yml      = createFileChannel(params.spades_yml)
    ch_spades_hmm      = createFileChannel(params.spades_hmm)
    ch_constraint_meta = createFileChannel(params.mapping_constraints)
    ch_annotation_meta = createFileChannel(params.annotation_metadata)

    // Databases, we really don't want to stage unnecessary databases
    ch_ref_pool      = (!params.skip_assembly && !params.skip_polishing) || (!params.skip_consensus_qc && !params.skip_blast_qc)           ? createChannel( params.reference_pool, "reference", true )                                                         : channel.empty()
    ch_kraken2_db    = (!params.skip_assembly && !params.skip_polishing && !params.skip_precluster) || !params.skip_read_classification    ? createChannel( params.kraken2_db, "kraken2", ('kraken2' in read_classifiers || 'kraken2' in contig_classifiers) ) : channel.empty()
    ch_kaiju_db      = (!params.skip_assembly && !params.skip_polishing && !params.skip_precluster) || !params.skip_read_classification    ? createChannel( params.kaiju_db, "kaiju", ('kaiju' in read_classifiers || 'kaiju' in contig_classifiers) )         : channel.empty()
    ch_checkv_db     = !params.skip_consensus_qc                                                                                           ? createChannel( params.checkv_db, "checkv", !params.skip_checkv )                                                  : channel.empty()
    ch_bracken_db    = !params.skip_read_classification                                                                                    ? createChannel( params.bracken_db, "bracken", ('bracken' in read_classifiers) )                                    : channel.empty()
    ch_k2_host       = !params.skip_preprocessing                                                                                          ? createChannel( params.host_k2_db, "k2_host", !params.skip_hostremoval )                                           : channel.empty()
    ch_annotation_db = !params.skip_consensus_qc                                                                                           ? createChannel( params.annotation_db, "annotation", !params.skip_consensus_annotation )                            : channel.empty()
    ch_prokka_db     = !params.skip_consensus_qc                                                                                           ? createChannel( params.prokka_db, "prokka", !params.skip_prokka )                                                  : channel.empty()

    // Importing samplesheet
    ch_reads = ch_samplesheet

    // Prepare Databases
    ch_db = channel.empty()
    if ((!params.skip_assembly && !params.skip_polishing) || !params.skip_consensus_qc || !params.skip_read_classification || (!params.skip_preprocessing && !params.skip_hostremoval)){

        ch_db_raw = ch_db.mix(ch_ref_pool,ch_kraken2_db, ch_kaiju_db, ch_checkv_db, ch_bracken_db, ch_k2_host, ch_annotation_db, ch_prokka_db)
        UNPACK_DB (ch_db_raw)

        ch_db = UNPACK_DB.out.db
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
                prokka: meta.id == 'prokka'
                    return [ unpacked ]
            }

        // transfer to value channels so processes are not just done once
        // '.collect()' is necessary to transform to list so cartesian products are made downstream
        ch_ref_pool         = ch_db.reference.collect{_meta, unpacked -> unpacked}.ifEmpty([]).map{unpacked -> [[id: 'reference'], unpacked]}
        ch_annotation_db    = ch_db.annotation.collect{_meta, unpacked -> unpacked}.ifEmpty([]).map{unpacked -> [[id: 'annotation'], unpacked]}
        ch_kraken2_db       = ch_db.kraken2.collect().ifEmpty([])
        ch_kaiju_db         = ch_db.kaiju.collect().ifEmpty([])
        ch_checkv_db        = ch_db.checkv.collect().ifEmpty([])
        ch_bracken_db       = ch_db.bracken.collect().ifEmpty([])
        ch_k2_host          = ch_db.k2_host.collect().ifEmpty([])
        ch_prokka_db        = ch_db.prokka.collect().ifEmpty([])
    }

    // Prepare blast DB
    ch_blast_refdb  = channel.empty()

    if ( params.reference_pool && ((!params.skip_assembly && !params.skip_polishing) || (!params.skip_consensus_qc && !params.skip_blast_qc))){
        BLAST_MAKEBLASTDB ( ch_ref_pool, [] )
        ch_blast_refdb = BLAST_MAKEBLASTDB.out.db
    }

    ch_host_trim_reads      = channel.empty()
    ch_decomplex_trim_reads = channel.empty()

    // preprocessing illumina reads
    if (!params.skip_preprocessing){
        PREPROCESSING_ILLUMINA (
            ch_reads,
            ch_k2_host,
            params.adapter_fasta ? file(params.adapter_fasta, checkIfExists:true) : [],
            ch_contaminants
            )
        ch_host_trim_reads      = PREPROCESSING_ILLUMINA.out.reads
        ch_decomplex_trim_reads = PREPROCESSING_ILLUMINA.out.reads_decomplexified
        ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.mqc.collect{_meta, mqc -> mqc}.ifEmpty([]))
        ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.low_reads_mqc.ifEmpty([]))
    } else {
        // Nothing downstream drops empty samples, so remove them here.
        // countFastq() reads every file on the head node, hence only when preprocessing is skipped.
        ch_host_trim_reads      = ch_reads.filter{ _meta, reads -> reads[0].countFastq() > 0}
        ch_decomplex_trim_reads = ch_reads.filter{ _meta, reads -> reads[0].countFastq() > 0}
    }

    // Reads used for downstream mapping & polishing steps (iterative refinement, mapping-constraint
    // selection, final variant-calling mapping).
    ch_mapping_polishing_reads = use_host_filtered_reads ? ch_host_trim_reads : ch_decomplex_trim_reads

    // Determining metagenomic diversity
    if (!params.skip_read_classification) {
        FASTQ_KRAKEN_KAIJU(
            ch_host_trim_reads,
            read_classifiers,
            ch_kraken2_db,
            ch_bracken_db,
            ch_kaiju_db,
            params.kraken2_save_reads,
            params.kraken2_save_readclassification,
            params.kaiju_taxon_rank
            )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_KRAKEN_KAIJU.out.mqc.collect{_meta, mqc -> mqc}.ifEmpty([]))
    }

    // Assembly
    ch_unaligned_raw_contigs     = channel.empty()
    ch_unaligned_contigs         = channel.empty()
    ch_polishing_consensus_reads = channel.empty()

    // channel for consensus sequences that have been generated across different iteration
    ch_consensus                 = channel.empty()
    // channel for consensus sequences that have been generated at the LAST iteration
    ch_consensus_reads           = channel.empty()
    // channel for summary table of clusters to include in mqc report
    ch_clusters_summary          = channel.empty()
    // channel for summary coverages of each contig
    ch_clusters_tsv              = channel.empty()

    if (!params.skip_assembly) {
        // run different assemblers and combine contigs
        FASTQ_ASSEMBLY( ch_host_trim_reads, ch_spades_yml, ch_spades_hmm, normalise_reads)
        ch_contigs       = FASTQ_ASSEMBLY.out.scaffolds
        ch_coverages     = FASTQ_ASSEMBLY.out.coverages
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ASSEMBLY.out.mqc.ifEmpty([]))

        if (!params.skip_polishing){
            // blast contigs against reference & identify clusters of (contigs & references)
            ch_contigs_reads = ch_contigs
                .join(ch_host_trim_reads, by: [0], remainder: false)

            FASTA_CONTIG_CLUST (
                ch_contigs_reads,
                ch_coverages,
                ch_blacklist,
                ch_blast_refdb.ifEmpty([]),
                ch_ref_pool.ifEmpty([]),
                ch_kraken2_db,
                ch_kaiju_db,
                contig_classifiers,
                params.cluster_method,
                params.identity_threshold,
                params.skip_precluster,
                params.perc_reads_contig,
                params.cluster_with_reference_pool
                )

            // Split up clusters into singletons and clusters of multiple contigs
            ch_centroids_members = FASTA_CONTIG_CLUST.out.centroids_members
                .map { meta, centroids, members ->
                    [ meta, centroids, members ]
                }
                .branch { meta, centroids, members ->
                    singletons: meta.cluster_size == 0
                        return [ meta + [step:"singleton"], centroids ]
                    multiple: meta.cluster_size >   0
                        return [ meta + [step:"consensus"], centroids, members ]
                }

            ch_clusters_summary    = FASTA_CONTIG_CLUST.out.clusters_summary.collect{_meta, summary -> summary}.ifEmpty([])
            ch_clusters_tsv        = FASTA_CONTIG_CLUST.out.clusters_tsv.collect{_meta, tsv -> tsv}.ifEmpty([])
            ch_multiqc_files       =  ch_multiqc_files.mix(FASTA_CONTIG_CLUST.out.no_blast_hits_mqc.ifEmpty([]))

            // map clustered contigs & create a single consensus per cluster
            ALIGN_COLLAPSE_CONTIGS (
                ch_centroids_members.multiple
                )

            SINGLETON_FILTERING (
                ch_centroids_members.singletons,
                params.min_contig_size,
                params.max_n_perc
                )

            ch_consensus = ALIGN_COLLAPSE_CONTIGS.out.consensus.mix( SINGLETON_FILTERING.out.filtered )

            ch_unaligned_raw_contigs = ALIGN_COLLAPSE_CONTIGS.out.unaligned_fasta
                .mix( SINGLETON_FILTERING.out.filtered )

            // We want the meta from the reference channel to be used downstream as this is our varying factor
            // To do this we combine the channels based on sample
            // Extract the reference meta's and reads
            // Make cartesian product of identified references & reads so all references will be mapped against again.
                ch_consensus_reads_intermediate = ch_consensus
                    .map { meta, fasta -> [meta.sample, meta, fasta] }
                    .combine(ch_mapping_polishing_reads.map { meta, fastq -> [meta.sample, meta, fastq]}, by: [0])
                    .map{
                        _sample, meta_ref, fasta, _meta_reads, fastq -> [meta_ref, fasta, fastq]
                    }

            if (!params.skip_iterative_refinement) {
                FASTQ_FASTA_ITERATIVE_CONSENSUS (
                    ch_consensus_reads_intermediate,
                    params.iterative_refinement_cycles,
                    params.intermediate_mapper,
                    params.with_umi,
                    params.deduplicate,
                    params.call_intermediate_variants,
                    params.intermediate_variant_caller,
                    params.intermediate_consensus_caller,
                    params.intermediate_mapping_stats,
                    params.min_mapped_reads,
                    keep_unmapped,
                    params.min_contig_size,
                    params.max_n_perc
                )
                ch_consensus                 = ch_consensus.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_allsteps)
                ch_polishing_consensus_reads = FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_reads
                ch_multiqc_files             = ch_multiqc_files.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.mqc.ifEmpty([])) //collect already done in subworkflow
            } else {
                ch_polishing_consensus_reads = ch_consensus_reads_intermediate
            }
        }
    }

    // add last step to it
    ch_consensus_reads = ch_polishing_consensus_reads
        .map{ meta, fasta, fastq ->
            [meta + [step: "variant-calling", iteration:'variant-calling', previous_step: meta.step], fasta, fastq]
        }

    ch_mash_screen = channel.empty()

    if (params.mapping_constraints && !params.skip_variant_calling ) {
        // Importing samplesheet
        ch_mapping_constraints = channel
            .fromList(samplesheetToList(params.mapping_constraints, "${projectDir}/assets/schemas/mapping_constraints.json"))
            .map { meta, sequence ->
                def samples = meta.samples == [] ? null : tuple(meta.samples.split(";"))  // Split up samples if meta.samples is not null
                [meta, samples, sequence]
            }
            .transpose(remainder: true)                                                   // Unnest

        // Joining all the reads with the mapping constraints, filter for those specified or keep everything if none specified.
        ch_map_seq_anno_combined = ch_mapping_polishing_reads
            .combine ( ch_mapping_constraints )
            .filter { meta_reads, _fastq, _meta_mapping, mapping_samples, _sequence ->
                mapping_samples == null || mapping_samples == meta_reads.sample
            }
            .map { meta, reads, meta_mapping, _samples, sequence_mapping ->
                def id = "${meta.sample}_${meta_mapping.id}-CONSTRAINT"
                def new_meta = meta + meta_mapping + [
                    id: id,
                    cluster_id: "${meta_mapping.id}",
                    step: "constraint",
                    isConstraint: true,
                    reads: reads,
                    iteration: 'variant-calling',
                    previous_step: 'constraint'
                    ]
                return [new_meta, sequence_mapping]
            }

        // Map with both reads and mapping constraints
        ch_constraint_consensus_reads = ch_map_seq_anno_combined
            .map { meta, fasta -> [meta, fasta, meta.reads] }
            .branch { meta, _fasta, _fastq ->
                multiFastaSelection : meta.selection == true
                singleFastaSelection : meta.selection == false
            }

        // Select the correct reference
        FASTQ_FASTA_MASH_SCREEN (
            ch_constraint_consensus_reads.multiFastaSelection
        )
        ch_mash_screen = FASTQ_FASTA_MASH_SCREEN.out.json.collect{_meta, json -> json}

        // For QC we keep original sequence to compare to
        ch_unaligned_contigs = ch_unaligned_raw_contigs
            .mix(ch_constraint_consensus_reads.singleFastaSelection.map { meta, fasta, _reads -> [meta, fasta] })
            .mix(FASTQ_FASTA_MASH_SCREEN.out.reference_fastq.map { meta, fasta, _reads -> [meta, fasta] })

        // Add to the consensus channel, which will be used for variant calling
        ch_consensus_reads = ch_consensus_reads
            .mix(FASTQ_FASTA_MASH_SCREEN.out.reference_fastq)
            .mix(ch_constraint_consensus_reads.singleFastaSelection)
    }

    // After consensus sequences have been made, we still have to map against it and call variants
    if ( !params.skip_variant_calling ) {

        FASTQ_FASTA_MAP_CONSENSUS(
            ch_consensus_reads,
            params.mapper,
            params.with_umi,
            params.deduplicate,
            true,
            params.variant_caller,
            params.consensus_caller,
            params.mapping_stats,
            params.min_mapped_reads,
            keep_unmapped,
            params.min_contig_size,
            params.max_n_perc
        )
        ch_consensus     = ch_consensus.mix(FASTQ_FASTA_MAP_CONSENSUS.out.consensus_all)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTA_MAP_CONSENSUS.out.mqc.ifEmpty([])) // collect already done in subworkflow

    }

    // Annotate variants with supplied gff
    if ( !params.skip_variant_calling && !params.skip_vcf_annotation ) {
        VCF_ANNOTATE (
            FASTQ_FASTA_MAP_CONSENSUS.out.vcf_ref
        )
    }

    ch_checkv_summary     = channel.empty()
    ch_quast_summary      = channel.empty()
    ch_blast_summary      = channel.empty()
    ch_annotation_summary = channel.empty()

    // Run several sumarisation tools on the variant calling results
    if ( !params.skip_consensus_qc  && (!params.skip_assembly || !params.skip_variant_calling) ) {
        ch_consensus_filter = ch_consensus
            .filter{_meta, fasta -> getLengthAndAmbigous(fasta).contig_size > 0}
        CONSENSUS_QC(
            ch_consensus_filter,
            ch_unaligned_contigs,
            ch_checkv_db,
            ch_blast_refdb.ifEmpty([]),
            ch_annotation_db,
            ch_prokka_db
            )
        ch_checkv_summary     = CONSENSUS_QC.out.checkv.collect{_meta, summary -> summary}.ifEmpty([])
        ch_quast_summary      = CONSENSUS_QC.out.quast.collect{_meta, summary -> summary}.ifEmpty([])
        ch_blast_summary      = CONSENSUS_QC.out.blast.collect{_meta, summary -> summary}.ifEmpty([])
        ch_annotation_summary = CONSENSUS_QC.out.annotation.collect{_meta, summary -> summary}.ifEmpty([])
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config           = params.multiqc_config              ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    // ch_multiqc_logo                    = params.multiqc_logo                ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_multiqc_custom_table_headers       = params.custom_table_headers        ? channel.fromPath(params.custom_table_headers, checkIfExists:true ) : channel.fromPath("$projectDir/assets/custom_table_headers.yml", checkIfExists:true )
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'viralmetagenome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // Prepare MULTIQC custom tables
    CUSTOM_MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_clusters_summary.ifEmpty([]),
        ch_metadata,
        ch_checkv_summary.ifEmpty([]),
        ch_quast_summary.ifEmpty([]),
        ch_blast_summary.ifEmpty([]),
        ch_constraint_meta,
        ch_annotation_summary.ifEmpty([]),
        ch_annotation_meta,
        ch_clusters_tsv.ifEmpty([]),
        ch_mash_screen.ifEmpty([]),
        ch_multiqc_custom_table_headers.ifEmpty([])
    )

    emit:
    multiqc_report = CUSTOM_MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
