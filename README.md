<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-viralgenie_logo_dark.png">
    <img alt="nf-core/viralgenie" src="docs/images/nf-core-viralgenie_logo_light.png">
  </picture>
</h1>

<!--[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/viralgenie/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/) [![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/Joon-Klaps/viralgenie)

[![GitHub Actions CI Status](https://github.com/Joon-Klaps/viralgenie/actions/workflows/ci.yml/badge.svg)](https://github.com/Joon-Klaps/viralgenie/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/Joon-Klaps/viralgenie/actions/workflows/linting.yml/badge.svg)](https://github.com/Joon-Klaps/viralgenie/actions?query=workflow%3A%22nf-core+linting%22)

<!-- [![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23viralgenie-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/viralgenie)-->

> [!TIP]
> Make sure to checkout the [viralgenie website](https://joon-klaps.github.io/viralgenie/) for more elaborate documentation!

## Introduction

**Viralgenie** is a bioinformatics best-practice analysis pipeline for reconstructing consensus genomes and to identify intra-host variants from metagenomic sequencing data or enriched based sequencing data like hybrid capture.

## Pipeline summary

![viralgenie-workflow](docs/images/metromap_style_pipeline_workflow_viralgenie.png)

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Performs optional read pre-processing
    - Adapter trimming([`fastp`](https://github.com/OpenGene/fastp), [`Trimmomatic`](https://github.com/usadellab/Trimmomatic))
    - Low complexity and quality filtering ([`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/))
    - Host-read removal ([`BowTie2`](http://bowtie-bio.sourceforge.net/bowtie2/))
3. Metagenomic diveristy mapping
    - Performs taxonomic classification and/or profiling using one or more of:
        - [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
        - [`Bracken`](https://ccb.jhu.edu/software/bracken/)[optional]
        - [`Kaiju`](https://kaiju.binf.ku.dk/)
    - Plotting Kraken2 and Kaiju ([`Krona`](https://hpc.nih.gov/apps/kronatools.html))
4. Denovo assembly ([`SPAdes`](http://cab.spbu.ru/software/spades/), [`TRINITY`](https://github.com/trinityrnaseq/trinityrnaseq), [`megahit`](https://github.com/voutcn/megahit)), combine contigs.
5. Contig reference idententification ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
    -   Identify top 5 blast hits
    -   Merge blast hit and all contigs of a sample
6. [Optional] Precluster contigs based on taxonomy
    - Identify taxonomy [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) and\or [`Kaiju`](https://kaiju.binf.ku.dk/)
    - Resolve potential inconsistencies in taxonomy & taxon filtering | simplification `bin/extract_precluster.py`
7. Cluster contigs (or every taxonomic bin) of samples, options are:
    - [`cdhitest`](https://sites.google.com/view/cd-hit)
    - [`vsearch`](https://github.com/torognes/vsearch/wiki/Clustering)
    - [`mmseqs-linclust`](https://github.com/soedinglab/MMseqs2/wiki#linear-time-clustering-using-mmseqs-linclust)
    - [`mmseqs-cluster`](https://github.com/soedinglab/MMseqs2/wiki#cascaded-clustering)
    - [`vRhyme`](https://github.com/AnantharamanLab/vRhyme)
    - [`Mash`](https://github.com/marbl/Mash)
8. Scaffolding of contigs to centroid ([`Minimap2`](https://github.com/lh3/minimap2), [`iVar-consensus`](https://andersen-lab.github.io/ivar/html/manualpage.html))
9. [Optional] Annotate 0-depth regions with external reference `bin/lowcov_to_reference.py`.
10. [Optional] Select best reference from `--mapping_constrains`:
    - [`Mash sketch`](https://github.com/marbl/Mash)
    - [`Mash screen`](https://github.com/marbl/Mash)
11. Mapping filtered reads to supercontig and mapping constrains([`BowTie2`](http://bowtie-bio.sourceforge.net/bowtie2/),[`BWAmem2`](https://github.com/bwa-mem2/bwa-mem2) and [`BWA`](https://github.com/lh3/bwa))
12. [Optional] Deduplicate reads ([`Picard`](https://broadinstitute.github.io/picard/) or if UMI's are used [`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html))
13. Variant calling and filtering ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html),[`iVar`](https://andersen-lab.github.io/ivar/html/manualpage.html))
14. Create consensus genome ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html),[`iVar`](https://andersen-lab.github.io/ivar/html/manualpage.html))
15. Repeat step 11-14 multiple times for the denovo contig route
16. Consensus evaluation and annotation ([`QUAST`](http://quast.sourceforge.net/quast),[`CheckV`](https://bitbucket.org/berkeleylab/checkv/src/master/),[`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [`mmseqs-search`](https://github.com/soedinglab/MMseqs2/wiki#batch-sequence-searching-using-mmseqs-search))
17. Result summary visualisation for raw read, alignment, assembly, variant calling and consensus calling results ([`MultiQC`](http://multiqc.info/))

## Usage

!!! Note
    If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Now, you can run the pipeline using:

```bash
nextflow run Joon-Klaps/viralgenie \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

!!! Warning
     Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
     see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://github.io/Joon-klaps/viralgenie/usage) and the [parameter documentation](https://github.io/Joon-klaps/viralgenie/parameters).

## Credits

Viralgenie was originally written by [`Joon-Klaps`](https://github.com/Joon-Klaps).

We thank the following people for their extensive assistance in the development of this pipeline:

-   [`Philippe Lemey`](https://github.com/plemey)
-   [`Liana Kafetzopoulou`](https://github.com/LianaKafetzopoulou)
-   [`nf-core community`](https://nf-co.re/)


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](https://github.io/Joon-klaps/viralgenie/CONTRIBUTING).

<!--
For further information or help, don't hesitate to get in touch on the [Slack `#viralgenie` channel](https://nfcore.slack.com/channels/viralgenie) (you can join with [this invite](https://nf-co.re/join/slack)).
-->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->

<!-- If you use  Joon-Klaps/viralgenie for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](https://github.io/Joon-klaps/viralgenie/CITATIONS) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
