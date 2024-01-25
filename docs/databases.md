# Joon-Klaps/viralgenie: Databases

## Introduction

Viralgenie uses a multitude of databases in order to analyse reads, contigs and consensus constructs. The default databases will be sufficient in most cases but there are always exceptions. This document will guide you towards the right documentation location for creating your custom databases.

> [!NOTE]
> Keep an eye out for [nf-core createtaxdb](https://nf-co.re/createtaxdb/) as it can be used for the main databases.

## Kaiju

<!-- TODO -->

## Reference pool

<!-- TODO -->

## Kraken2 databases

A number of database indexes have already been generated and maintained by [@BenLangmead Lab](https://github.com/BenLangmead), see [here](https://benlangmead.github.io/aws-indexes/k2). These databases can directly be used to run the workflow with Kraken2 as well as Bracken.

In case the databases above do not contain your desired libraries, you can build a custom Kraken2 database. This requires two components: a taxonomy (consisting of `names.dmp`, `nodes.dmp`, and `*accession2taxid`) files, and the FASTA files you wish to include.

To pull the NCBI taxonomy, you can run the following:

```bash
kraken2-build --download-taxonomy --db <YOUR_DB_NAME>
```

You can then add your FASTA files with the following build command.

```bash
kraken2-build --add-to-library *.fna --db <YOUR_DB_NAME>
```

You can repeat this step multiple times to iteratively add more genomes prior building.

Once all genomes are added to the library, you can build the database (and optionally clean it up):

```bash
kraken2-build --build --db <YOUR_DB_NAME>
kraken2-build --clean --db <YOUR_DB_NAME>
```

You can then add the `<YOUR_DB_NAME>/` path to your nf-core/taxprofiler database input sheet.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `kraken2`
  - `opts.k2d`
  - `hash.k2d`
  - `taxo.k2d`

</details>

You can follow the Kraken2 [tutorial](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) for a more detailed description.

### Host read removal

Viralgenie uses kraken2 to remove contaminated reads.

> motivation read: [The human “contaminome”: bacterial, viral, and computational contamination in whole genome sequences from 1000 families](https://www.nature.com/articles/s41598-022-13269-z) and [Reconstruction of the personal information from human genome reads in gut metagenome sequencing data](https://www.nature.com/articles/s41564-023-01381-3)

The contamination database is likely the largest database. The default databases is made small explicitly made smaller to save storage for end users but is not optimal. I would recommend to create a database consisting of the libraries `human, archea, bacteria` which will be more then 200GB in size. Additionally, it's good practice to include DNA & RNA of the host of origin if known (i.e. mice, ticks, mosquito, ... ). Add it as described above.

Set it with the variable `--host_k2_db`

### Viral Diversity with Kraken2

The metagenomic diveristy estimated with kraken2 is based on the viral refseq database which can cut short if you expect your the species within your sample to have a large amount of diversity eg below 80% ANI ([quasi-species](https://link.springer.com/chapter/10.1007/978-3-642-77011-1_1)). To resolve this it's better to create a database that contains a wider species diversity then only one genome per species. Databases that have this wider diversity is [Virosaurus](https://viralzone.expasy.org/8676) or the [RVDB](https://rvdb.dbi.udel.edu/home) which can increase the accuracy of kraken2. Add it as described above.

Set ut with the variable kraken2_db `--kraken2_db`

### Bracken
