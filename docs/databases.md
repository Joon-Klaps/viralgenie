# Joon-Klaps/viralgenie: Databases

## Introduction

Viralgenie uses a multitude of databases in order to analyse reads, contigs and consensus constructs. The default databases will be sufficient in most cases but there are always exceptions. This document will guide you towards the right documentation location for creating your custom databases.

> [!NOTE]
> Keep an eye out for [nf-core createtaxdb](https://nf-co.re/createtaxdb/) as it can be used for the main databases.

## Reference pool

The reference pool dataset is used to identify potential references for scaffolding. It's a fasta file that will be used to make a blast database within the pipeline. The default database is the [clustered Reference Viral DataBase (C-RVDB)](https://rvdb.dbi.udel.edu/) a database that was build for enhancing virus detection using high-throughput/next-generation sequencing (HTS/NGS) technologies. An alternative reference pool is the [Virosaurus](https://viralzone.expasy.org/8676) database which is a manually curated database of viral genomes.

Any nucleotide fasta file will do.

Specify it with the parameter `--reference_pool`.

## Annotation sequences

Identifying the species and the segment of the final genome construct is done based on a blast search to a smaller annotated sequencing data `assets/bv-brc-refvirus-anno.fasta.gz`. This dataset is extracted from the [BV-BRC](https://www.bv-brc.org/). The default annotation dataset was reconstructed using their [CLI-tool](https://www.bv-brc.org/docs/cli_tutorial/index.html) to communicate with their [relational database](https://www.bv-brc.org/docs/cli_tutorial/cli_getting_started.html#the-bv-brc-database)

In case you want to modify the default database see [here how to install the CLI](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html) and use the following scripts as a guideline to reconstruct the dataset.

Here we select all viral genomes that are not lab reassortments and are reference genomes and add metadata attributes to the output.

```bash
# download annotation metadata +/- 5s
p3-all-genomes --eq superkingdom,Viruses --eq reference_genome,Reference --ne host_common_name,'Lab reassortment' --attr genome_id,species,segment,genome_name,genome_length,host_common_name,genbank_accessions,taxon_id   > all-virus-anno.txt
# download genome data, done seperatly as it takes much longer to query +/- 1 hour
p3-all-genomes --eq superkingdom,Viruses --eq reference_genome,Reference --ne host_common_name,'Lab reassortment' | p3-get-genome-contigs --attr sequence > all-virus.fasta
```

> [!NOTE]
> Any attribute can be downloaded and will be added to the final report if the formatting remains the same.
> For a complete list of attributes see `p3-all-genomes --fields` or read their [manual](https://www.bv-brc.org/docs/cli_tutorial/cli_getting_started.html)

Next, the metadata and the genomic data is combined into a single fasta file where the metada fields are stored in the fasta comment as `key1:"value1"|key2:"value2"|...` using the following python code.

```python
import pandas as pd
import re

# read in sequences with all columns as strings
sequences = pd.read_csv("refseq-virus.fasta", index_col=0, sep="\t", dtype=str)
data = pd.read_csv("refseq-virus-anno.txt", index_col=0, sep="\t", dtype=str)

# merge the df's
df = sequences.join(data)
# remove 'genome' from the name
df.columns = df.columns.str.replace("genome.", "")

# create fasta header
def create_fasta_header(row):
    annotations = "|".join(
        [
            f'{column}:"{value}"'
            for column, value in row.items()
            if column != "contig.sequence"
        ]
    )
    return f"{annotations}\n"


df["fasta_header"] = df.apply(create_fasta_header, axis=1)

df["fasta_entry"] = (
    ">" + df.index.astype(str) + " " + df["fasta_header"] + df["contig.sequence"]
)
with open("bv-brc-refvirus-anno.fasta", "w") as f:
    for entry in df["fasta_entry"]:
        f.write(entry + "\n")
```

<details markdown="1">
<summary>Expected files in database directory</summary>

-   `bv-brc-refvirus-anno.fasta.gz`

</details>

This annotation database can then be specified using `--annotation_db`

## Kaiju

The kaiju database will be used to classify the reads and intermediate contigs in taxonomic groups. The default database is the RVDB-prot pre-built database from Kaiju.

A number of Kaiju pre-built indexes for reference datasets are maintained by the the developers of Kaiju and made available on the [Kaiju website](https://bioinformatics-centre.github.io/kaiju/downloads.html).
To build a kaiju database, you need three components: a FASTA file with the protein sequences, the NCBI taxonomy dump files, and you need to define the uppercase characters of the standard 20 amino acids you wish to include.

> [!Warning]
> The headers of the protein fasta file must be numeric NCBI taxon identifiers of the protein sequences.

To download the NCBI taxonomy files, please run the following commands:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```

To build the database, run the following command (the contents of taxdump must be in the same location where you run the command):

```bash
kaiju-mkbwt -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
kaiju-mkfmi proteins
```

> [!Tip]
> You can speed up database construction by supplying the threads parameter (`-t`).

<details markdown="1">
<summary>Expected files in database directory</summary>

-   `kaiju`
    -   `kaiju_db_*.fmi`
    -   `nodes.dmp`
    -   `names.dmp`

</details>

For the Kaiju database construction documentation, see [here](https://github.com/bioinformatics-centre/kaiju#custom-database).

## Kraken2 databases

The Kraken2 database will be used to classify the reads and intermediate contigs in taxonomic groups.

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

-   `kraken2`
    -   `opts.k2d`
    -   `hash.k2d`
    -   `taxo.k2d`

</details>

You can follow the Kraken2 [tutorial](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) for a more detailed description.

### Host read removal

Viralgenie uses kraken2 to remove contaminated reads.

> motivation read: [The human “contaminome”: bacterial, viral, and computational contamination in whole genome sequences from 1000 families](https://www.nature.com/articles/s41598-022-13269-z) and [Reconstruction of the personal information from human genome reads in gut metagenome sequencing data](https://www.nature.com/articles/s41564-023-01381-3)

The contamination database is likely the largest database. The default databases is made small explicitly made smaller to save storage for end users but is not optimal. I would recommend to create a database consisting of the libraries `human, archea, bacteria` which will be more then 200GB in size. Additionally, it's good practice to include DNA & RNA of the host of origin if known (i.e. mice, ticks, mosquito, ... ). Add it as described above.

Set it with the variable `--host_k2_db`

### Viral Diversity with Kraken2

The metagenomic diveristy estimated with kraken2 is based on the viral refseq database which can cut short if you expect your the species within your sample to have a large amount of diversity eg below 80% ANI ([quasi-species](https://link.springer.com/chapter/10.1007/978-3-642-77011-1_1)). To resolve this it's better to create a database that contains a wider species diversity then only one genome per species. Databases that have this wider diversity is [Virosaurus](https://viralzone.expasy.org/8676) or the [RVDB](https://rvdb.dbi.udel.edu/home) which can increase the accuracy of kraken2. Add it as described above.

Set it with the variable `--kraken2_db`

### Bracken

Bracken does not require an independent database but rather builds upon Kraken2 databases. [The pre-built Kraken2 databases hosted by Ben Langmead](https://benlangmead.github.io/aws-indexes/k2) already contain the required files to run Bracken.

However, to build custom databases, you will need a Kraken2 database, the (average) read lengths (in bp) of your sequencing experiment, the K-mer size used to build the Kraken2 database, and Kraken2 available on your machine.

```bash
bracken-build -d <KRAKEN_DB_DIR> -k <KRAKEN_DB_KMER_LENGTH> -l <READLENGTH>
```

> [!Tip]
> You can speed up database construction by supplying the threads parameter (`-t`).

> [!Tip]
> If you do not have Kraken2 in your `$PATH` you can point to the binary with `-x /<path>/<to>/kraken2`.

<details markdown="1">
<summary>Expected files in database directory</summary>

-   `bracken`
    -   `hash.k2d`
    -   `opts.k2d`
    -   `taxo.k2d`
    -   `database100mers.kmer_distrib`
    -   `database150mers.kmer_distrib`

</details>

You can follow Bracken [tutorial](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual) for more information.

Set it with the variable `--bracken_db`
