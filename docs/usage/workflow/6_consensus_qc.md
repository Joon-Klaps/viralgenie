# Report generation and quality control

nf-core/viralmetagenome's report and result interpretation heavily relies on MultiQC. MultiQC is a tool to create a single report from multiple analysis results. It is designed to be used with a wide range of bioinformatics tools and is compatible with a wide range of data formats. Almost all tools are summarised within the MultiQC report that have interactive plots and data tables. However, due to the number of tools included, some results are summarised in the directory `overview-tables` to reduce the size of the MultiQC report.

:::tip
Complete output descriptions of files and images can be found in the [output section](../output.md).
:::

Within the MultiQC report, nf-core/viralmetagenome provides a number of custom tables based on consensus genome quality control data. These tools are:

- [QUAST](#quast): QUAST is a quality assessment tool for genome assemblies. It calculates various metrics such as N50, L50, number of contigs, and total length.
- [CheckV](#checkv): CheckV is a tool for assessing the quality of metagenome-assembled viral genomes. It calculates various metrics such as completeness, contamination, and strain heterogeneity.
- [Prokka](#prokka): Prokka is a whole genome annotation pipeline for identifying features of interest in a set of genomic DNA sequences, and labelling them with useful information.
- [blastn](#blast): BLAST is a tool for comparing primary biological sequence information. It calculates the similarity between the consensus genome and the reference genome.
- [mmseqs-search](#mmseqs-search) - included as 'annotation': MMseqs is an ultra-fast and sensitive search tool for protein and nucleotide databases. nf-core/viralmetagenome uses MMseqs to annotate the consensus genomes and assign them a species name, segment name, expected host, etc.
- [mafft](#mafft): MAFFT is a multiple sequence alignment program.
- [SnpEff and SnpSift](#snpeff-and-snpsift): SnpEff is a genetic variant annotation and functional effect prediction tool. SnpSift is a toolbox that allows you to filter and manipulate annotated files.

> [!NOTE]
> Consensus genome quality control can be skipped with `--skip_consensus_qc`.

## QUAST

[QUAST](http://quast.sourceforge.net/quast) is a quality assessment tool for genome assemblies. It calculates various metrics such as N50, L50, number of contigs, and total length. However, in the summary table, it is mainly used to get the number of ambiguous bases in the consensus genome.

> [!NOTE]
> QUAST can be skipped with `--skip_quast`.

## CheckV

[CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) is a tool for assessing the quality of metagenome-assembled viral genomes. It calculates various metrics such as completeness, contamination, and strain heterogeneity. CheckV estimates completeness by comparing sequences with a large database of complete viral genomes, metagenomes, metatranscriptomes, and metaviromes.

:::tip{title="Incomplete genomes for segmented viruses"}
CheckV estimates the completeness of a virus based on all genome segments. If a virus has multiple segments, the completeness of the virus is calculated based on the length of the concatenated segments. For example, Lassa virus has 2 segments L: 7.2kb and S: 3.4kb. The completeness of the virus is calculated based on the length of the concatenated segments (7.2kb + 3.4kb = 10.6kb) and so if the generated consensus genome of the L segment is 7.1kb it will report the completeness as 7.1/10.6 ~ 67%.
:::

> [!NOTE]
> CheckV can be skipped with `--skip_checkv`.

## Prokka

[Prokka](https://github.com/tseemann/prokka) is a whole genome annotation pipeline for identifying features of interest in a set of genomic DNA sequences, and labelling them with useful information. Prokka is a software tool to annotate bacterial, archaeal and viral genomes.

:::tip{title="Suboptimal annotation"}
Prokka was initially designed for bacterial and archaeal genomes, and may not be optimal for viral genomes. [VIGOR4](https://github.com/JCVenterInstitute/VIGOR4) is a good alternative but is species specific.
:::

:::tip{title="Custom protein database"}
Prokka can be given a custom protein database to annotate your genomes with, have a look at [prot-RVDB](https://rvdb-prot.pasteur.fr/) for viral protein databases. Supply the database using `--prokka_db`.
:::

> [!NOTE]
> Prokka can be skipped with `--skip_prokka`.

## BLAST

[blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is a tool for comparing primary biological sequence information. It calculates the similarity between the consensus genome and the reference genome. The similarity is calculated based on the number of identical bases between the two sequences.nf-core/viralmetagenome uses blastn to compare the sequences against the supplied `--reference_pool` dataset.

> [!NOTE]
> BLASTn can be skipped with `--skip_blast_qc`.

## MMseqs-search

[MMseqs-search](https://github.com/soedinglab/MMseqs2/wiki#searching) is an ultra-fast and sensitive search tool for protein and nucleotide databases.nf-core/viralmetagenome uses MMseqs to search the consensus genomes in an annotated database, like [Virosaurus](https://virosaurus.vital-it.ch/) (see also [defining your own custom annotation database](../databases.md#annotation-sequences)), and uses the annotation data of the best hit to assign the consensus genome a species name, segment name, expected host, and any other metadata that describes the database sequences - taken from a metadata table when `--annotation_metadata` is given, and from the fasta headers otherwise. This allows nf-core/viralmetagenome, in addition to the BLAST search of reference pool hits, to compare the generated consensus genomes at a species & segment level.

:::info
MMseqs was used for the annotation step instead of BLAST because of the ability to query using a tblastx search for highly diverging viruses while supplying a nucleotide annotation database. To specify another type of search (e.g. blastp, blastx, etc.), please refer to the [parameters consensus-qc section](../parameters.md#consensus-qc).
:::

> [!NOTE]
> MMseqs-search can be skipped with `--skip_consensus_annotation`.

## SnpEff and SnpSift

[SnpEff](https://pcingola.github.io/SnpEff/) is a genetic variant annotation and functional effect prediction tool. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

[SnpSift](https://pcingola.github.io/SnpEff/SnpSift.html) is a toolbox that allows you to filter and manipulate annotated files. The ExtractFields tool is used to extract specific information from the annotated VCF files into a tabular format for easier analysis.

nf-core/viralmetagenome uses SnpEff to annotate variants identified by the variant calling process with functional information, and SnpSift ExtractFields to extract key information from the annotated variants into a more accessible tabular format.

The annotation process provides valuable information about the impact of variants, including:

- Whether variants are synonymous or non-synonymous
- Changes in amino acid sequences
- Potential impact severity (HIGH, MODERATE, LOW, MODIFIER)
- Gene and transcript information

> [!NOTE]
> Variant annotation can be skipped with `--skip_vcf_annotation`.

## MAFFT

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is a multiple sequence alignment program. It is used to align the following genomic data:

- The final consensus genome
- The identified reference genome from `--reference_pool`
- The de novo contigs from each assembler (that constituted the final consensus genome)
- Each consensus genome from the iterative refinement steps.

> [!NOTE]
> MAFFT can be skipped with `--skip_alignment_qc`.

## MultiQC

[MultiQC](https://multiqc.info/) is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.

<image src="https://raw.githubusercontent.com/MultiQC/MultiQC/main/docs/images/multiqc_overview.excalidraw.svg"/>
> Image credit: [MultiQC](https://multiqc.info/)

Reports are generated by scanning given directories for recognised log files. These are parsed and a single HTML report is generated summarising the statistics for all logs found. MultiQC reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control.

MultiQC is also used to generate the `overview-tables` as it extracts the additional data from various tools. The data that needs to be extracted can be modified with the argument `--custom_table_headers` where a [yml file](https://github.com/nf-core/viralmetagenome/blob/dev/assets/custom_table_headers.yml) shows which tools need to be included in the summary table in addition to BLAST, CheckV, QUAST, and MMseqs (annotation).

```yml title="custom_table_headers.yml"
tool:
  - Tool subsection: # if applicable
      - name_in_mqc_table: "new name"
      - output_reads: "deduplicated reads"
      - percent_passing_dedup: "% passing dedup"
```
