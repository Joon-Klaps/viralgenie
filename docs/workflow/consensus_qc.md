# Report generation and quality control

Viralgenie's report and result interpretation heavily relies on MultiQC. MultiQC is a tool to create a single report from multiple analysis results. It is designed to be used with a wide range of bioinformatics tools and is compatible with a wide range of data formats. Almost all tools are summarised within the MultiQC report that have interactive plots and data tables. However, due to the number of tools included, some results are summarised in the directory `overview-tables` to reduce the size of the MultiQC report.

!!! Tip
    Complete output descriptions of files and images can be found in the [output section](../output.md).

Within the MultiQC report, Viralgenie provides a number of custom tables based on consensus genome quality control data. These tools are:

- [QUAST](#quast): QUAST is a quality assessment tool for genome assemblies. It calculates various metrics such as N50, L50, number of contigs, and total length.
- [CheckV](#checkv): CheckV is a tool for assessing the quality of metagenome-assembled viral genomes. It calculates various metrics such as completeness, contamination, and strain heterogeneity.
- [blastn](#blastn): BLAST is a tool for comparing primary biological sequence information. It calculates the similarity between the consensus genome and the reference genome.
- [mmseqs-search](#mmseqs-search) - included as 'annotation': MMseqs is an ultra-fast and sensitive search tool for protein and nucleotide databases. Viralgenie uses MMseqs to annotate the consensus genomes and assign them a species name, segment name, expected host, etc.
- [mafft](#mafft): MAFFT is a multiple sequence alignment program.

> Consensus genome quality control can be skipped with `--skip_consensus_qc`.

## QUAST

[QUAST](http://quast.sourceforge.net/quast) is a quality assessment tool for genome assemblies. It calculates various metrics such as N50, L50, number of contigs, and total length. However, in the summary table, it is mainly used to get the number of ambiguous bases in the consensus genome.

> QUAST can be skipped with `--skip_quast`.

## CheckV

[CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) is a tool for assessing the quality of metagenome-assembled viral genomes. It calculates various metrics such as completeness, contamination, and strain heterogeneity. CheckV estimates completeness by comparing sequences with a large database of complete viral genomes, metagenomes, metatranscriptomes, and metaviromes.

!!! Tip "Incomplete genomes for segmented viruses"
    CheckV estimates the completeness of a virus based on all genome segments. If a virus has multiple segments, the completeness of the virus is calculated based on the length of the concatenated segments. For example, Lassa virus has 2 segments L: 7.2kb and S: 3.4kb. The completeness of the virus is calculated based on the length of the concatenated segments (7.2kb + 3.4kb = 10.6kb) and so if the generated consensus genome of the L segment is 7.1kb it will report the completeness as 7.1/10.6 ~ 67%.

> CheckV can be skipped with `--skip_checkv`.

## BLASTn

[blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is a tool for comparing primary biological sequence information. It calculates the similarity between the consensus genome and the reference genome. The similarity is calculated based on the number of identical bases between the two sequences. Viralgenie uses blastn to compare the sequences against the supplied `--reference_pool` dataset.

> BLASTn can be skipped with `--skip_blast_qc`.

## MMseqs-search

[MMseqs-search](https://github.com/soedinglab/MMseqs2/wiki#searching) is an ultra-fast and sensitive search tool for protein and nucleotide databases. Viralgenie uses MMseqs to search the consensus genomes in an annotated database, like [Virousarus](https://virosaurus.vital-it.ch/) (see also [defining your own custom annotation database](../customisation/databases.md#annotation-sequences)), and uses the annotation data of the best hit to assign the consensus genome a species name, segment name, expected host, and any other metadata that is embedded within the database. This allows Viralgenie, in addition to the BLAST search of reference pool hits, to compare the generated consensus genomes at a species & segment level.

!!! info
    MMseqs was used for the annotation step instead of BLAST because of the ability to query using a tblastx search for highly diverging viruses while supplying a nucleotide annotation database. To specify another type of search (e.g. blastp, blastx, etc.), please refer to the [parameters consensus-qc section](../parameters.md#consensus-qc).

> MMseqs-search can be skipped with `--skip_annotation`.


## MAFFT

[MAFFT](https://mafft.cbrc.jp/alignment/software/) is a multiple sequence alignment program. It is used to align the following genomic data:

- The final consensus genome
- The identified reference genome from `--reference_pool`
- The de novo contigs from each assembler (that constituted the final consensus genome)
- Each consensus genome from the iterative refinement steps.

> MAFFT can be skipped with `--skip_alignment_qc`.

## MultiQC

[MultiQC](https://multiqc.info/) is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.

<image src="https://raw.githubusercontent.com/MultiQC/MultiQC/main/docs/images/multiqc_overview.excalidraw.svg"/>

Reports are generated by scanning given directories for recognised log files. These are parsed and a single HTML report is generated summarising the statistics for all logs found. MultiQC reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control.

MultiQC is also used to generate the `overview-tables` as it extracts the additional data from various tools. The data that needs to be extracted can be modified with the argument `--custom_table_headers` where a [yml file](https://github.com/Joon-Klaps/viralgenie/blob/dev/assets/custom_table_headers.yml) shows which tools need to be included in the summary table in addition to BLAST, CheckV, QUAST, and MMseqs (annotation).

```yml title="custom_table_headers.yml"
tool:
  - Tool subsection : # if applicable
    - name_in_mqc_table: "new name"
    - output_reads: "deduplicated reads"
    - percent_passing_dedup: "% passing dedup"
```
