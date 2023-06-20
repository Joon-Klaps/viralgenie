# ![nf-core/viralgenie](docs/images/nf-core-viralgenie_logo_light.png#gh-light-mode-only) ![nf-core/viralgenie](docs/images/nf-core-viralgenie_logo_dark.png#gh-dark-mode-only)

<!--[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/viralgenie/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<!--
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/viralgenie)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23viralgenie-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/viralgenie)
-->

[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/viralgenie** is a bioinformatics best-practice analysis pipeline for reconstructing consensus denovo genomes and identify intra-host variants from metagenomic sequencing data or enriched based sequencing data like hybrid capture.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/viralgenie/results).

## Pipeline summary

**Pipeline is still under active development and not all steps have been implemented.**

![viral-genie-workflow](docs/images/workflow-v3.png)

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
<!-- TODO nf-core: include the [`???`] in all steps-->

Tools given for each step can change as benchmarking is in progress.

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),[`falco`](https://github.com/smithlabcode/falco))
2. Performs optional read pre-processing
   - Adapter trimming([`fastp`](https://github.com/OpenGene/fastp), [`Trimmomatic`](https://github.com/usadellab/Trimmomatic))
   - Low complexity and quality filtering ([`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/))
   - Host-read removal ([`BowTie2`](http://bowtie-bio.sourceforge.net/bowtie2/))
3. Supports statistics for host-read removal ([`Samtools`](http://www.htslib.org/))
4. Metagenomic diveristy mapping
   - Performs taxonomic classification and/or profiling using one or more of:
     - [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
     - [`Kaiju`](https://kaiju.binf.ku.dk/)
   - Perform optional post-processing ([`bracken`](https://ccb.jhu.edu/software/bracken/))
   - Plotting Kraken2, Centrifuge, Kaiju and MALT results ([`Krona`](https://hpc.nih.gov/apps/kronatools.html))
5. Denovo assembly ([`SPAdes`](http://cab.spbu.ru/software/spades/), [`TRINITY`](https://github.com/trinityrnaseq/trinityrnaseq), [`megahit`](https://github.com/voutcn/megahit)), combine contigs.
6. Taxonomy annotation and contig reference identification ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))

- Identify top 5 blast hits
- Merge blast hit seqs and contigs into a single file & cluster ([`cdhit`](https://sites.google.com/view/cd-hit), [`vsearch`](https://github.com/torognes/vsearch/wiki/Clustering))
- Filter for only contigs that correspond to the identified centroid

7. Scaffolding of contigs ([`AGA`](https://github.com/emweb/aga), [`minimap2`](https://github.com/lh3/minimap2), [`MUMMER`](https://github.com/mummer4/mummer)
8. Creation of hybrid supercontig by combining contigs and the reference.
9. Mapping filtered reads to hybrid supercontig ([`BowTie2`](http://bowtie-bio.sourceforge.net/bowtie2/),[`minimap2`](https://github.com/lh3/minimap2)
10. Update consensus genome with mapped reads ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html),[`GATK`](https://github.com/broadinstitute/gatk))
11. Repeat step 9-10 multiple times
12. Variant calling and annotation ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html),[`GATK`](https://github.com/broadinstitute/gatk),[`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
13. Assembly & contig evaluation ([`QUAST`](http://quast.sourceforge.net/quast),[`CheckV`](https://bitbucket.org/berkeleylab/checkv/src/master/))
14. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/viralgenie \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://nf-co.re/viralgenie/usage) and the [parameter documentation](https://nf-co.re/viralgenie/parameters).

## Pipeline output

The nf-core/viralgenie pipeline comes with documentation about the pipeline [usage](https://nf-co.re/viralgenie/usage), [parameters](https://nf-co.re/viralgenie/parameters) and [output](https://nf-co.re/viralgenie/output).
-->

## Credits

viralgenie was originally written by joon.klaps@kuleuven.be.

We thank the following people for their extensive assistance in the development of this pipeline:

- [`Philippe Lemey`](https://github.com/plemey)
- [`Liana Kafetzopoulou`](https://github.com/LianaKafetzopoulou)

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

<!--
For further information or help, don't hesitate to get in touch on the [Slack `#viralgenie` channel](https://nfcore.slack.com/channels/viralgenie) (you can join with [this invite](https://nf-co.re/join/slack)).
-->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/viralgenie for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
