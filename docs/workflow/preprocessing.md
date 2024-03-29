# Preprocessing

Viralgenie offers three main preprocessing steps for the preprocessing of raw sequencing reads:

- [Read quality control](#read-quality-control): read quality assessment and filtering.
- [Read processing](#read-processing): adapter clipping and pair-merging.
- [Complexity filtering](#complexity-filtering): removal of low-sequence complexity reads.
- [Host read-removal](#host-read-removal): removal of reads aligning to reference genome(s) of a host.

![preprocessing](../images/preprocessing.png)

> Preprocessing can be entirely skipped with the option `--skip_preprocessing`.
> See the [parameters preprocessing section](../parameters.md#preprocessing-options) for all relevant arguments to control the preprocessing steps.

## Read Quality control
[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used before and after read processing and after host read-removal to assess the quality of the reads.

```mermaid
graph LR;
    A[Raw reads] --> B["`**FastQC**`"];
    B --> C[Read processing];
    C --> D["`**FastQC**`"];
    D --> E[Complexity filtering];
    E --> G[Host read-removal];
    G --> H["`**FastQC**`"];
```

## Read processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging is performed by the tools [`fastp`](https://github.com/OpenGene/fastp) or [`Trimmomatic`](https://github.com/usadellab/Trimmomatic). The tool `fastp` is a fast all-in-one tool for preprocessing fastq files. The tool `Trimmomatic` is a flexible read trimming tool for Illumina NGS data. Both tools can be used to remove adapters and low-quality reads from the raw sequencing reads. An adapter file can be provided through the argument `--adapter_fasta`.

Specify the tool to use for read processing with the `--trim_tool` parameter, the default is `fastp`.

## Complexity filtering

Complexity filtering is primarily a run-time optimisation step. Low-complexity sequences are defined as having commonly found stretches of nucleotides with limited information content (e.g. the dinucleotide repeat CACACACACA). Such sequences can produce a large number of high-scoring but biologically insignificant results in database searches. Removing these reads therefore saves computational time and resources.

Complexity filtering is done with [`Bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) which is part of [`BBtools`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) where the "duk" stands for Decontamination Using Kmers.

By default this step is skipped, if this step shouldn't be skipped specify `--skip_complexity_filtering false`.

## Host read-removal

To avoid the inclusion of contamined reads into consensus genomes or when determining intra-host variants, host read-removal is performed. Host read-removal is performed by the tool `Kraken2`.

!!! info
    The reason why we use Kraken2 for host removal over regular read mappers is nicely explained in the following papers:

    * [The human “contaminome”: bacterial, viral, and computational contamination in whole genome sequences from 1000 families](https://www.nature.com/articles/s41598-022-13269-z)
    * [Reconstruction of the personal information from human genome reads in gut metagenome sequencing data](https://www.nature.com/articles/s41564-023-01381-3)

Specify the host database with the `--host_k2_db` parameter. The default is small subset of the human genome and we highly suggest that you make this database more elaborate. For this, read the section on [creating custom kraken2 host databases](../customisation/databases.md#kraken2-databases).
