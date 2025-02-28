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
| `with_umi` | With or without UMI detection |  |
| `skip_umi_extract` | With or without UMI extraction | True |
| `umi_deduplicate` | Specify at what level UMI deduplication should occur. | read |
| `umi_discard_read` | Discard R1 / R2 if required 0, meaning not to discard | 0 |
| `trim_tool` | The used trimming tool | fastp |
| `skip_trimming` | Skip read trimming |  |
| `adapter_fasta` | Fasta file of adapters |  |
| `save_trimmed_fail` | Specify true to save files that failed to pass trimming thresholds ending in `*.fail.fastq.gz` |  |
| `save_merged` | Specify true to save all merged reads to a file ending in `*.merged.fastq.gz` |  |
| `min_trimmed_reads` | Inputs with fewer than this reads will be filtered out of the "reads" output channel | 1 |
| `skip_complexity_filtering` | Skip filtering of low complexity regions in reads <details><summary>Help</summary><small>Low-complexity sequences are defined as having commonly found stretches of nucleotides with limited information content (e.g. the dinucleotide repeat CACACACACA). Such sequences can produce a large number of high-scoring but biologically insignificant results in database searches</small></details>| True |
| `decomplexifier` | Specify the decomplexifier to use, bbduk or prinseq | prinseq |
| `contaminants` | Reference files containing adapter and/or contaminant sequences for sequence kmer matching (used by bbduk) |  |
| `skip_hostremoval` | Skip the removal of host read sequences |  |
| `host_k2_db` | Kraken2 database used to remove host and contamination | s3://ngi-igenomes/test-data/viralrecon/kraken2_human.tar.gz |
| `host_k2_library` | Kraken2 library(s) required to remove host and contamination <details><summary>Help</summary><small>Only used when no host kraken2 database is specified.</small></details>| human |
| `skip_host_fastqc` | Skip the fastqc step after host & contaminants were removed |  |
| `arguments_fastqc` | Arguments for FastQC tool | --quiet |
| `arguments_fastp` | Arguments for Fastp tool | --cut_front --cut_tail --trim_poly_x --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 10 --length_required 50 |
| `arguments_trimmomatic` | Arguments for Trimmomatic tool | ILLUMINACLIP:null:2:30:10 |
| `arguments_umitools_extract` | Arguments for UMI-tools extract |  |
| `arguments_humid` | Arguments for Humid tool | -a -m 1 |
| `arguments_bbduk` | Arguments for BBDuk tool | entropy=0.3 entropywindow=50 entropymask=f |
| `arguments_prinseq_reads` | Arguments for Prinseq tool for reads |  |
| `arguments_kraken2_host` | Arguments for Kraken2 tool for host removal |  |

## Metagenomic diversity

Parameters used to determine the metagenomic diversity of the sample

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_read_classification` | Skip determining the metagenomic diversity of the sample |  |
| `read_classifiers` | Specify the taxonomic read classifiers, choices are 'kaiju,kraken2' | kraken2,kaiju |
| `save_databases` | Save the used databases |  |
| `kraken2_db` | Location of the Kraken2 database | https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230314.tar.gz |
| `kraken2_save_reads` | Save classified and unclassified reads as fastq files |  |
| `kraken2_save_readclassification` | Save summary overview of read classifications in a txt file |  |
| `kraken2_save_minimizers` | Save kraken2's used minimizers |  |
| `bracken_db` | Location of bracken database | https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20230314.tar.gz |
| `kaiju_db` | Location of Kaiju database | https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_rvdb_2023-05-26.tgz |
| `kaiju_taxon_rank` | Level of taxa rank that needs to be determined | species |
| `arguments_kraken2` | Arguments for Kraken2 tool | --report-minimizer-data |
| `arguments_kaiju` | Arguments for Kaiju tool | -v |
| `arguments_kaiju2table` | Arguments for Kaiju2Table tool | -e -l species |
| `arguments_kaiju2krona` | Arguments for Kaiju2Krona tool | -v -u |
| `arguments_krona` | Arguments for Krona tool |  |
| `arguments_bracken` | Arguments for Bracken tool |  |
| `arguments_kreport2krona` | Arguments for Kreport2Krona tool |  |

## Assembly

Parameters relating to the used assembly methods

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_assembly` | Skip de novo assembly of reads |  |
| `assemblers` | The specified tools for de novo assembly, multiple options are possible | spades,megahit |
| `spades_mode` | Specific SPAdes mode to run | rnaviral |
| `spades_hmm` | File or directory with amino acid HMMs for Spades HMM-guided mode. |  |
| `spades_yml` | Path to yml file containing read information. <details><summary>Help</summary><small>The raw FASTQ files listed in this YAML file MUST be supplied to the respective illumina/pacbio/nanopore input channel(s) _in addition_ to this YML. File entries in this yml must contain only the file name and no paths.</small></details>|  |
| `skip_contig_prinseq` | Skip the filtering of low complexity contigs with prinseq |  |
| `skip_sspace_basic` | Skip the contig extension with sspace_basic |  |
| `read_distance` | Specify the mean distance between the paired reads | 350 |
| `read_distance_sd` | Specify the deviation of the mean distance that is allowed. <details><summary>Help</summary><small>For instance, a mean of 200 and a sd of 0.75. This means that any pair having a distance between 150 and 250 is allowed.</small></details>| 0.75 |
| `read_orientation` | Specify the read orientation. | FR |
| `arguments_spades` | Arguments for SPAdes tool | --rnaviral |
| `arguments_megahit` | Arguments for MEGAHIT tool |  |
| `arguments_trinity` | Arguments for Trinity tool | --max_reads_per_graph 100000 |
| `arguments_quast` | Arguments for QUAST tool | --min-contig 0 |
| `arguments_sspace_basic` | Arguments for SSPACE Basic tool | -x 1 -o 15 -r 0.75 |
| `arguments_prinseq_contig` | Arguments for Prinseq tool for contigs | -out_format 1 -lc_dust .20 |

## Polishing

Parameters relating to the refinement of de novo contigs

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_polishing` | Skip the refinement/polishing of contigs through reference based scaffolding and read mapping |  |
| `save_intermediate_polishing` | Save intermediate polishing files <details><summary>Help</summary><small>There are multiple processes within the polishing subworkflow that might not contain relevant information</small></details>|  |
| `reference_pool` | Set of fasta sequences used as potential references for the contigs | https://rvdb.dbi.udel.edu/download/C-RVDBvCurrent.fasta.gz |
| `skip_precluster` | Skip the preclustering of assemblies to facilitate downstream processing of assemblies |  |
| `keep_unclassified` | Keep the contigs that could not be classified with the taxonomic databases (`kaiju_db` & `kraken2_db`) <details><summary>Help</summary><small>Within the preclustering step, all contigs will get a taxonomic classification using the provided databases for the metagenomic tools. In some cases, the number of unclassified contigs can be very large if the database is restrictive. This will result in large clusters in downstream processing that can take up a lot of resources despite not being a priority in some analyses. So set it to `True` if you want to keep unclassified contigs and set it to `False` if you don't want to keep them. </small></details>| True |
| `precluster_classifiers` | Specify the metagenomic classifiers to use for contig taxonomy classification: 'kraken2,kaiju' | kraken2,kaiju |
| `cluster_method` | Cluster algorithm used for contigs | cdhitest |
| `network_clustering` | (only with mash) Algorithm to partition the network. <details><summary>Help</summary><small>Mash creates a distance matrix that gets translated into a network of connectected nodes where the edges represent the similarity. This network is then split up using the specified method.<br><br> - [leiden](https://leidenalg.readthedocs.io/en/stable/intro.html) algorithm: a hierarchical clustering algorithm, that recursively merges communities into single nodes by greedily optimizing the modularity<br> - [connected_components] algorithm: a clustering algorithm that defines the largest possible communities where each node within a subset is reachable from every other node in the same subset via any edge .<br></small></details>| connected_components |
| `skip_nocov_to_reference` | Skip creation of the hybrid consensus, instead keep the scaffold with ambiguous bases if the depth of scaffolds is not high enough. |  |
| `identity_threshold` | Identity threshold value used in clustering algorithms | 0.85 |
| `perc_reads_contig` | Minimum cumulated sum of mapped read percentages of each member from a cluster group, set to 0 to disable <details><summary>Help</summary><small>Setting this variable will remove clusters that have a low cumulated sum of mapped read percentages. This can be used to remove clusters that have a low coverage and are likely to be false positives.</small></details>| 5 |
| `min_contig_size` | Minimum allowed contig size <details><summary>Help</summary><small>Setting this to a low value will result in a large number of questionable contigs and an increase in computation time </small></details>| 500 |
| `max_contig_size` | Maximum allowed contig size | 10000000 |
| `max_n_perc` | Define the maximum percentage of ambiguous bases in a contig | 50 |
| `skip_singleton_filtering` | Skip the filtering of contigs that did not cluster together with other contigs <details><summary>Help</summary><small>Setting this to true will cause the pipeline not to remove contigs that don't have similar contigs. Filtering settings can be further specified with `min_contig_size` and `max_n_100kbp`.</small></details>|  |
| `arguments_blast_makeblastdb` | Arguments for BLAST makeblastdb tool | -dbtype nucl |
| `arguments_blastn` | Arguments for BLASTN tool | -max_target_seqs 5 |
| `arguments_blast_filter` | Arguments for BLAST filter tool | --escore 0.01 --bitscore 50 --percent-alignment 0.80 |
| `arguments_kraken2_contig` | Arguments for Kraken2 tool for contigs |  |
| `arguments_kaiju_contig` | Arguments for Kaiju tool for contigs | -v |
| `arguments_extract_precluster` | Arguments for precluster extraction | --keep-unclassified true --merge-strategy lca |
| `arguments_cdhit` | Arguments for CD-HIT tool | -c 0.85 -mask rRyYkKsSwWmMbBdDhHvVnN |
| `arguments_vsearch` | Arguments for VSEARCH tool | --maxseqlength 10000000 --id 0.85 --strand both --iddef 0 --no_progress --qmask none |
| `arguments_mmseqs_linclust` | Arguments for MMseqs2 linclust tool | --min-seq-id 0.85 -c 0.700 --cov-mode 2 --cluster-mode 0 |
| `arguments_mmseqs_cluster` | Arguments for MMseqs2 cluster tool | --min-seq-id 0.85 -c 0.700 --cov-mode 2 --cluster-mode 0 |
| `arguments_vrhyme` | Arguments for VRhyme tool | --mems 50 |
| `arguments_mash_dist` | Arguments for Mash distance tool | -s 4000 -k 15 |
| `arguments_network_cluster` | Arguments for network clustering | --score 0.85 |
| `arguments_extract_cluster` | Arguments for cluster extraction | --perc_reads_contig 5 |
| `arguments_minimap2_align` | Arguments for Minimap2 alignment |  |
| `arguments_minimap2_index` | Arguments for Minimap2 index |  |
| `arguments_mash_sketch` | Arguments for Mash sketch tool | -i |
| `arguments_mash_screen` | Arguments for Mash screen tool |  |
| `arguments_select_reference` | Arguments for selecting reference |  |

## Iterative consensus refinement

Define parameters for iterations to update de novo consensus using  reference based improvements

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_iterative_refinement` | Don't realign reads to consensus sequences and redefine the consensus through (multiple) iterations |  |
| `iterative_refinement_cycles` | Number of iterations | 2 |
| `intermediate_mapper` | Mapping tool used during iterations | bwamem2 |
| `intermediate_variant_caller` | Variant caller used during iterations | ivar |
| `call_intermediate_variants` | Call variants during the iterations <details><summary>Help</summary><small>Will always be done when iterative consensus caller is bcftools</small></details>|  |
| `intermediate_consensus_caller` | Consensus tool used for calling new consensus during iterations | bcftools |
| `intermediate_mapping_stats` | Calculate summary statistics during iterations | True |

## Variant analysis

Parameters relating to the analysis of variants associated to contigs and scaffolds

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_variant_calling` | Skip the analysis of variants for the external reference or contigs |  |
| `mapper` | Define which mapping tool needs to be used when mapping reads to reference | bwamem2 |
| `mapping_constraints` | Sequence to use as a mapping reference instead of the de novo contigs or scaffolds |  |
| `deduplicate` | Deduplicate the reads <details><summary>Help</summary><small>If used with UMI's, `umi tools` will be used to group and call consensus of each individual read group. If not used with UMI's use `PicardsMarkDuplicates`. </small></details>| True |
| `variant_caller` | Define the variant caller to use: 'ivar' or 'bcftools' | ivar |
| `consensus_caller` | Consensus tool used for calling new consensus in final iteration | ivar |
| `min_mapped_reads` | Define the minimum number of mapped reads in order to continue the variant and consensus calling | 200 |
| `allele_frequency` | Minimum allele frequency threshold for calling consensus | 0.75 |
| `mapping_stats` | Calculate summary statistics in final iteration | True |
| `ivar_header` |  |  |
| `arguments_bwamem2_index` | Arguments for BWA-MEM2 index |  |
| `arguments_bwa_index` | Arguments for BWA index |  |
| `arguments_bwa_mem` | Arguments for BWA MEM |  |
| `arguments_bowtie2_build` | Arguments for Bowtie2 build |  |
| `arguments_bowtie2_align` | Arguments for Bowtie2 alignment | --local --very-sensitive-local --seed 1 |
| `arguments_umitools_dedup` | Arguments for UMI-tools deduplication | --umi-separator=\':\' --method cluster --unmapped-reads use |
| `arguments_picard_markduplicates` | Arguments for Picard MarkDuplicates | --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp --REMOVE_DUPLICATES true |
| `arguments_picard_collectmultiplemetrics` | Arguments for Picard CollectMultipleMetrics | --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp |
| `arguments_custom_mpileup` | Arguments for custom mpileup |  |
| `arguments_mosdepth` | Arguments for Mosdepth tool |  |
| `arguments_bcftools_mpileup1` | Arguments for BCFtools mpileup step 1 | --ignore-overlaps --count-orphans --max-depth 800000 --min-BQ 20 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR |
| `arguments_bcftools_mpileup2` | Arguments for BCFtools mpileup step 2 | --ploidy 2 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only |
| `arguments_bcftools_mpileup3` | Arguments for BCFtools mpileup step 3 | --include \'INFO/DP>=10\ |
| `arguments_bcftools_norm` | Arguments for BCFtools norm | --do-not-normalize --output-type z --multiallelics -any |
| `arguments_bcftools_stats` | Arguments for BCFtools stats |  |
| `arguments_samtools_stats` | Arguments for Samtools stats command |  |
| `arguments_samtools_idxstats` | Arguments for Samtools idxstats command |  |
| `arguments_samtools_flagstat` | Arguments for Samtools flagstat command |  |
| `arguments_tabix` | Arguments for Tabix tool | -p vcf -f |
| `arguments_bedtools_merge` | Arguments for Bedtools merge |  |
| `arguments_bedtools_maskfasta` | Arguments for Bedtools maskfasta |  |
| `arguments_bcftools_consensus` | Arguments for BCFtools consensus |  |
| `arguments_ivar_variants1` | Arguments for iVar variants step 1 | -q 20 -m 10 |
| `arguments_ivar_variants2` | Arguments for iVar variants step 2 | --ignore-overlaps --count-orphans --max-depth 0 --no-BAQ --min-BQ 0 |
| `arguments_make_bed_mask` | Arguments for making BED mask | -a --ignore-overlaps --count-orphans --max-depth 0 --no-BAQ --min-BQ 0 |
| `arguments_ivar_consensus1` | Arguments for iVar consensus step 1 | -t 0 -q 20 -m 10 -n N |
| `arguments_ivar_consensus2` | Arguments for iVar consensus step 2 | --count-orphans --max-depth 0 --min-BQ 20 --no-BAQ -aa |

## Consensus QC

Apply different quality control techniques on the generated consensus genomes

| Parameter | Description | Default |
|-----------|-----------|-----------|
| `skip_consensus_qc` | Skip the quality measurements on consensus genomes |  |
| `skip_checkv` | Skip the use of checkv for quality check |  |
| `checkv_db` | Reference database used by checkv for consensus quality control <details><summary>Help</summary><small>If not given, the most recent one is downloaded.</small></details>|  |
| `skip_annotation` | Skip the annotation of the consensus constructs |  |
| `annotation_db` | Database used for annotation of the consensus constructs <details><summary>Help</summary><small>The metadata fields are stored in the fasta comment as `key1:"value1"|key2:"value2"|...` see docs/databases.md for more information.</small></details>| ftp://ftp.expasy.org/databases/viralzone/2020_4/virosaurus90_vertebrate-20200330.fas.gz |
| `skip_prokka` | Skip gene estimation & annotation with prokka |  |
| `prokka_db` | Define a prokka `--protein` database for protein annotation <details><summary>Help</summary><small>Specify a custom protein database for Prokka annotation</small></details>|  |
| `skip_quast` | Skip the use of QUAST for quality check |  |
| `skip_blast_qc` | Skip the blast search of contigs to the provided reference DB |  |
| `skip_alignment_qc` | Skip creating an alignment of each the collapsed clusters and each iterative step | True |
| `mmseqs_searchtype` | Specify the search algorithm to use for mmseqs. 0: auto 1: amino acid, 2: translated, 3: nucleotide, 4: translated nucleotide alignment <details><summary>Help</summary><small>Only search-type 3 supports both forward and reverse search<br>1 - BLASTP;<br>2 - TBLASTN;<br>3 - BLASTN;<br>4 - TBLASTX</small></details>| 4 |
| `arguments_checkv` | Arguments for CheckV tool | --remove_tmp |
| `arguments_mafft_iterations` | Arguments for MAFFT iterations | --auto --adjustdirection |
| `arguments_mafft_qc` | Arguments for MAFFT QC | --auto --adjustdirection |
| `arguments_blastn_qc` | Arguments for BLASTN QC | -max_target_seqs 5 |
| `arguments_prokka` | Arguments for Prokka tool | --centre X --compliant --force --kingdom Viruses |
| `arguments_mmseqs_search` | Arguments for MMseqs2 search | --search-type 4 --rescore-mode 3 |
| `arguments_quast_qc` | Arguments for QUAST quality control |  |

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
| `custom_table_headers` | Custom yaml file containing the table column names selection and new names. | https://github.com/Joon-Klaps/viralgenie/raw/refs/heads/dev/assets/custom_table_headers.yml |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | True |
| `prefix` | Prefix of all output files followed by [date]_[pipelineversion]_[runName] <details><summary>Help</summary><small>Use '--global_prefix' to not have metadata embedded.</small></details>|  |
| `global_prefix` | Global prefix set if you don't want metadata embedded in the output filenames |  |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | https://raw.githubusercontent.com/nf-core/test-datasets/ |
