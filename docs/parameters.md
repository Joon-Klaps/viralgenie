---
hide:
  - navigation
---

# Viralgenie pipeline parameters

A pipeline to reconstruct consensus genomes and identify intrahost variants from metagenomic sequencing data or enriched based sequencing data like hybrid capture.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row. See [usage docs](https://joon-klaps.github.io/viralgenie/latest/usage/#input).</small></details>|  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. |  |
| `metadata` | Sample metadata that is included in the multiqc report |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>|  |

## Preprocessing options

Options related to the trimming, low complexity and host removal steps of the reads

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_preprocessing` | Skip read preprocessing and use input reads for downstream analysis |  |
| `skip_fastqc` | Skip read quality statistics summary tool 'fastqc' |  |
| `save_final_reads` | Save reads after the final preprocessing step | True |
| `save_intermediate_reads` | Save reads after every preprocessing step |  |
| `with_umi` | With or without umi detection |  |
| `skip_umi_extract` | With or without umi extraction | True |
| `umi_deduplicate` | Specify at what level UMI deduplication should occur. | read |
| `humid_mismatches` | Specify the maximum number of mismatches between reads for them to still be considered neighbors. | 1 |
| `humid_strategy` | Specify the strategy for umi-deduplication directional vs cluster | directional |
| `umitools_dedup_strategy` | Specify the strategy or method for umi-tools deduplication on mapping level | cluster |
| `umi_discard_read` | Discard R1 / R2 if required 0, meaning not to discard | 0 |
| `trim_tool` | The used trimming tool | fastp |
| `skip_trimming` | Skip read trimming |  |
| `fastp_deduplicate` | Use Fastp's deduplicate option |  |
| `fastp_dedup_accuracy` | Define the accuracy used for hashes while deduplicating with faspt |  |
| `adapter_fasta` | Fasta file of adapters |  |
| `save_trimmed_fail` | Specify true to save files that failed to pass trimming thresholds ending in `*.fail.fastq.gz` |  |
| `save_merged` | Specify true to save all merged reads to the a file ending in `*.merged.fastq.gz` |  |
| `min_trimmed_reads` | Inputs with fewer than this reads will be filtered out of the "reads" output channel | 1 |
| `skip_complexity_filtering` | Skip filtering of low complexity regions in reads <details><summary>Help</summary><small>Low-complexity sequences are defined as having commonly found stretches of nucleotides with limited information content (e.g. the dinucleotide repeat CACACACACA). Such sequences can produce a large number of high-scoring but biologically insignificant results in database searches</small></details>| True |
| `decomplexifier` | Specify the decomplexifier to use, bbduk or prinseq | prinseq |
| `contaminants` | Reference files containing adapter and/or contaminant sequences for sequence kmer matching (used by bbduk) |  |
| `skip_hostremoval` | Skip the removal of host read sequences |  |
| `host_k2_db` | Kraken2 database used to remove host and conamination | s3://ngi-igenomes/test-data/viralrecon/kraken2_human.tar.gz |
| `host_k2_library` | Kraken2 library(s) required to remove host and contamination <details><summary>Help</summary><small>Only used when no host kraken2 database is specified.</small></details>| human |
| `skip_host_fastqc` | Skip the fastqc step after host & contaminants were removed |  |

## Metagenomic diveristy

Parameters used to determine the metagenomic diversity of the sample

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_read_classification` | Skip determining the metagenomic diversity of the sample |  |
| `read_classifiers` | Specify the taxonomic read classifiers, choices are 'kaiju,kraken2' | kraken2,kaiju |
| `save_databases` | Save the used databases |  |
| `kraken2_db` | Location of the Kraken2 database | https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230314.tar.gz |
| `kraken2_save_reads` | Save classified and unclassified reads  as fastq files |  |
| `kraken2_save_readclassification` | Save summary overview of read classifications in a txt file |  |
| `kraken2_save_minimizers` | Save kraken2's used minimizers |  |
| `bracken_db` | Location of bracken database | https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230314.tar.gz |
| `kaiju_db` | Location of Kaiju database | https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_rvdb_2023-05-26.tgz |
| `kaiju_taxon_rank` | Level of taxa rank that needs to be determined | species |

## Assembly

Parameters relating to the used assembly methods

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_assembly` | Skip de novo assembly of reads |  |
| `assemblers` | The specified tools for denovo assembly, multiple options are possible | spades,megahit |
| `spades_mode` | specific SPAdes mode to run | rnaviral |
| `spades_hmm` | File or directory with amino acid HMMs for Spades HMM-guided mode. |  |
| `spades_yml` | Path to yml file containing read information. <details><summary>Help</summary><small>The raw FASTQ files listed in this YAML file MUST be supplied to the respective illumina/pacbio/nanopore input channel(s) _in addition_ to this YML. File entries in this yml must contain only the file name and no paths.</small></details>|  |
| `assembler_patterns` | Regex pattern to identify contigs that have been made by the assemblers |  |
| `skip_contig_prinseq` | Skip the filtering of low complexity contigs with prinseq |  |
| `skip_sspace_basic` | Skip the contig extension with sspace_basic |  |
| `read_distance` | Specify the mean distance between the paired reads | 350 |
| `read_distance_sd` | Specify the a deviation of the mean distance that is allowed. <details><summary>Help</summary><small>For instance, a mean of 200 and a sd of 0.75. This means that any pair having a distance between 150 and 250 is allowed.</small></details>| 0.75 |
| `read_orientation` | Specify the read orientation. | FR |

## Polishing

Parameters relating to the refinement of denovo contigs

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_polishing` | Skip the refinement/polishing of contigs through reference based scaffolding and read mapping |  |
| `save_intermediate_polishing` | Save intermediate polishing files <details><summary>Help</summary><small>There are multiple processes within the polishing subworkflow that might not contain relevant information  </small></details>|  |
| `reference_pool` | Set of fasta sequences used as potential references for the contigs | https://rvdb.dbi.udel.edu/download/C-RVDBvCurrent.fasta.gz |
| `skip_precluster` | Skip the preclustering of assemblies to facilitate downstream processing of assemblies |  |
| `keep_unclassified` | Keep the contigs that could not be classified with the taxonomic databases (`kaiju_db` & `kraken2_db`) <details><summary>Help</summary><small>Within the preclustering step, all contigs will get a taxonomic classification using the provided databases for the metagenomic tools. In some cases, the number of unclassified contigs, can be very large if the database is restrictive. This will result in large clusters in downstream processing that can take up a lot of resources despite not being a priority in some analyses. So set it to `True` if you want to keep unclassified contigs and set it to `False` if you don't want to keep them. </small></details>| True |
| `precluster_classifiers` | Specify the metagenomic classifiers to use for contig taxonomy classification: 'kraken2,kaiju' | kraken2,kaiju |
| `precluster_merge_strategy` | Taxon conflict resolution mode, must be 1 (Kaiju), 2 (Kraken),  lca, or lowest. <details><summary>Help</summary><small>The option -c determines the method of resolving conflicts in the taxonomic assignment for a read.<br>Possible values are '1', '2', 'lca', 'lowest':<br>  '1' -> the taxon id from Kaiju is used.<br>  '2' -> the taxon id from Kraken is used.<br>  'lca' -> the least common ancestor of the two taxon ids from both input files is used.<br>  'lowest' -> the lower rank of the two taxa is used if they are within the same lineage. Otherwise the LCA is used.</small></details>| lca |
| `precluster_simplification_level` | Level of taxonomic simplification <details><summary>Help</summary><small>The taxonomic classification of the contigs can be simplified to a certain level. This can be done by specifying the taxonomic rank to which the classification should be simplified.</small></details>|  |
| `precluster_exclude_taxa` | Hard constrain for taxa to exclude from the preclustering, if multiple given make sure to enclose with '"' and separate with a space. <details><summary>Help</summary><small>The taxonomic classification of the contigs can be filtered. This can be done by specifying a list of taxa that should be excluded from the results.</small></details>|  |
| `precluster_exclude_children` | Taxon ids to exclude along with all their children from the preclustering, if multiple given make sure to enclose with '"' and separate with a space. <details><summary>Help</summary><small>The taxonomic classification of the contigs can be filtered. This can be done by specifying a list of taxa that should be excluded along with all their children from the results.</small></details>|  |
| `precluster_exclude_parents` | Taxon ids to exclude along with all their parents from the preclustering, if multiple given make sure to enclose with '"' and separate with a space. <details><summary>Help</summary><small>The taxonomic classification of the contigs can be filtered. This can be done by specifying a list of taxa that should be excluded along with all their parents from the results.</small></details>|  |
| `precluster_include_children` | Taxon ids to include along with all their children from the preclustering, if multiple given make sure to enclose with '"' and separate with a space. <details><summary>Help</summary><small>The taxonomic classification of the contigs can be filtered. This can be done by specifying a list of taxa that should be included along with all their children from the results.</small></details>|  |
| `precluster_include_parents` | Taxon ids to include along with all their parents from the preclustering, if multiple given make sure to enclose with '"' and separate with a space. <details><summary>Help</summary><small>The taxonomic classification of the contigs can be filtered. This can be done by specifying a list of taxa that should be included along with all their parents from the results.</small></details>|  |
| `cluster_method` | Cluster algorithm used for contigs | cdhitest |
| `mmseqs_cluster_mode` | Specify the algorithm to partition the network graph from mmseqs <details><summary>Help</summary><small>The Greedy Set cover (0) algorithm is an approximation for the NP-complete optimization problem called set cover.<br>Connected component (1) uses transitive connection to cover more remote homologs.<br>Greedy incremental (2) works analogous to CD-HIT clustering algorithm.</small></details>| 0 |
| `network_clustering` | (only with mash) Algorithm to partition the network. <details><summary>Help</summary><small>Mash creates a distance matrix that gets translated into a network of connectected nodes where the edges represent the similarity. This network is then split up using the specified method.<br><br> - [leiden](https://leidenalg.readthedocs.io/en/stable/intro.html) algorithm: a hierarchical clustering algorithm, that recursively merges communities into single nodes by greedily optimizing the modularity<br> - [connected_components] algorithm: a clustering algorithm that defines the largest possible communities where each node within a subset is reachable from every other node in the same subset via any edge .<br></small></details>| connected_components |
| `skip_hybrid_consensus` | Skip creation of the hybrid consensus, instead keep the scaffold with ambiguous bases if the depth of scaffolds is not high enough. |  |
| `identity_threshold` | Identity threshold value used in clustering algorithms | 0.85 |
| `min_contig_size` | Minimum allowed contig size <details><summary>Help</summary><small>Setting this to a low value will result in a large number of questionable contigs and an increase in computation time </small></details>| 500 |
| `perc_reads_contig` | minimum cumulated sum of mapped read percentages of each member from a cluster group, set to 0 to disable <details><summary>Help</summary><small>Setting this variable will remove clusters that have a low cumulated sum of mapped read percentages. This can be used to remove clusters that have a low coverage and are likely to be false positives.</small></details>| 5 |
| `max_contig_size` | Maximum allowed contig size | 10000000 |
| `max_n_perc` | Define the maximum percentage of ambiguous bases in a contig | 50 |
| `skip_singleton_filtering` | Skip the filtering of contigs that did not cluster together with other contigs <details><summary>Help</summary><small>Setting this to true will cause the pipeline not to remove contigs that don't have similar contigs. Filtering settings can be further specified with `min_contig_size` and `max_n_100kbp`.</small></details>|  |

## Iterative consensus refinement

Define parameters for iterations to update denovo consensus using  reference based improvements

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_iterative_refinement` | Don't realign reads to consensus sequences and redefine the consensus through (multiple) iterations |  |
| `iterative_refinement_cycles` | number of iterations | 2 |
| `intermediate_mapper` | mapping tool used during iterations | bwamem2 |
| `intermediate_variant_caller` | variant caller used during iterations | ivar |
| `call_intermediate_variants` | call variants during the iterations <details><summary>Help</summary><small>Will always be done when iterative consensus caller is bcftools</small></details>|  |
| `intermediate_consensus_caller` | consensus tool used for calling new consensus during iterations | bcftools |
| `intermediate_mapping_stats` | calculate summary statistics during iterations | True |

## Variant analysis

Parameters relating to the analysis of variants associated to contigs and scaffolds

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_variant_calling` | Skip the analysis of variants for the external reference or contigs |  |
| `mapper` | Define which mapping tool needs to be used when mapping reads to reference | bwamem2 |
| `mapping_constraints` | Sequence to use as a mapping reference instead of the de novo contigs or scaffolds |  |
| `deduplicate` | deduplicate the reads <details><summary>Help</summary><small>If used with umi's, `umi tools` will be used to group and call consensus of each indiual read group. If not used with umi's use `PicardsMarkDuplicates`. </small></details>| True |
| `variant_caller` | Define the variant caller to use: 'ivar' or 'bcftools' | ivar |
| `consensus_caller` | consensus tool used for calling new consensus in final iteration | ivar |
| `umi_separator` | UMI seperator in fastq header. <details><summary>Help</summary><small>If you have used an alternative method which does not separate the read id and UMI with a “_”, such as bcl2fastq which uses “:”, you can specify the separator with the option --umi_separator=<sep>, replacing <sep> with e.g “:”.</small></details>| : |
| `mash_sketch_size` | Specify the sketch size, the number of (non-redundant) min-hashes that are kept. <details><summary>Help</summary><small>Larger sketches will better represent the sequence, but at the cost of larger sketch files and longer comparison times.</small></details>| 4000 |
| `mash_sketch_kmer_size` | Specify the kmer size for mash to create their hashes <details><summary>Help</summary><small>larger k-mers will provide more specificity, while smaller k-mers will provide more sensitivity. Larger genomes will also require larger k-mers to avoid k-mers that are shared by chance</small></details>| 15 |
| `min_mapped_reads` | Define the minimum number of mapped reads in order to continue the variant and consensus calling | 200 |
| `mapping_stats` | calculate summary statistics in final iteration | True |
| `ivar_header` |  |  |

## Consensus QC

Apply different quality control techniques on the generated consensus genomes

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_consensus_qc` | Skip the quality measurements on consensus genomes |  |
| `skip_checkv` | Skip the use of checkv for quality check |  |
| `checkv_db` | Reference database used by checkv for consensus quality control <details><summary>Help</summary><small>If not given, the most recent one is downloaded.</small></details>|  |
| `skip_annotation` | Skip the annotation of the consensus constructs |  |
| `annotation_db` | Database used for annotation of the cosensus constructs <details><summary>Help</summary><small>The metada fields are stored in the fasta comment as `key1:"value1"|key2:"value2"|...` see docs/databases.md for more information.</small></details>| ftp://ftp.expasy.org/databases/viralzone/2020_4/virosaurus90_vertebrate-20200330.fas.gz |
| `skip_quast` | Skip the use of QUAST for quality check |  |
| `skip_blast_qc` | Skip the blast search of contigs to the provided reference DB |  |
| `skip_alignment_qc` | Skip creating an alignment of each the collapsed clusters and each iterative step |  |
| `mmseqs_searchtype` | Specify the search algorithm to use for mmseqs. 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment <details><summary>Help</summary><small>Only search-type 3 supports both forward and reverse search<br>1 - BLASTP;<br>2 - TBLASTN;<br>3 - BLASTN;<br>4 - TBLASTX</small></details>| 4 |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | master |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| https://raw.githubusercontent.com/nf-core/configs/master |
| `config_profile_name` | Institutional config name. |  |
| `config_profile_description` | Institutional config description. |  |
| `config_profile_contact` | Institutional config contact information. |  |
| `config_profile_url` | Institutional config URL link. |  |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `version` | Display version and exit. |  |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| copy |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>|  |
| `plaintext_email` | Send plain-text email instead of HTML. |  |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | 25.MB |
| `monochrome_logs` | Do not use coloured log outputs. |  |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>|  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. |  |
| `multiqc_config` | Custom config file to supply to MultiQC. |  |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | https://github.com/Joon-Klaps/viralgenie/blob/dev/docs/images/ViralGenie-nf-core-theme.png?raw=true |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. |  |
| `clean_output_on_error` | Delete the output directory if the pipeline fails |  |
| `custom_table_headers` | Custom yaml file contian g the table column names selection and new names. | ${projectDir}/assets/custom_table_headers.yml |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | True |
| `prefix` | Prefix of all output files followed by [date]_[pipelineversion]_[runName] <details><summary>Help</summary><small>Use '--global_prefix' to not have metadata embedded.</small></details>|  |
| `global_prefix` | Global prefix set if you don't want metadata embedded in the prefix |  |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | https://raw.githubusercontent.com/nf-core/test-datasets/ |
