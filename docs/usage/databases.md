# Databases

## Introduction

nf-core/viralmetagenome uses a multitude of databases in order to analyze reads, contigs, and consensus constructs. The default databases will be sufficient in most cases but there are always exceptions. This document will guide you towards the right documentation location for creating your custom databases.

:::tip
Building custom databases for taxonomic profilers, can be challenging. [nf-core createtaxdb](https://nf-co.re/createtaxdb/) adresses this issue!
:::

## Reference pool

The reference pool dataset is used to identify potential references for scaffolding. It's a multifasta file of diverse viral genomes in nucleotide format, for which a blast database will be made within the pipeline.

The default database is the [v31.0 of the Reference Viral DataBase (C-RVDB; Jan 9, 2026)](https://rvdb.dbi.udel.edu/) a database that was built for enhancing virus detection using high-throughput/next-generation sequencing (HTS/NGS) technologies. The RVDB is updated biannually.

An alternative reference pool is the [Virosaurus](https://viralzone.expasy.org/8676) database which is a manually curated database of viral genomes, which was last updated in 2020.

> [!NOTE]
> Some reference pools bundle partial or defective genomes (for example specific entries in RVDB). If you encounter such sequences, supply the `--blacklist` parameter with a newline-separated list of identifiers or identifier fragments to keep them out of the scaffolding step.

Any nucleotide fasta file will do. Specify it with the parameter `--reference_pool`.

> [!NOTE]
> The input to `--reference_pool` must be a multifasta file (can be compressed `.gz`), directories or `tar.gz` will fail.

## Kaiju

The Kaiju database will be used to classify the reads and intermediate contigs in taxonomic groups. The default database is the RVDB-prot pre-built database from Kaiju.

A number of Kaiju pre-built indexes for reference datasets are maintained by the developers of Kaiju and made available on the [Kaiju website](https://bioinformatics-centre.github.io/kaiju/downloads.html).
To build a Kaiju database, you need three components: a FASTA file with the protein sequences, the NCBI taxonomy dump files, and you need to define the uppercase characters of the standard 20 amino acids you wish to include.

:::warning
The headers of the protein fasta file must be numeric NCBI taxon identifiers of the protein sequences.

For example, a valid FASTA header would look like this:

```text
>1_1358
MAQQRRGGFKRRKKVDFIAANKIEVVDYKDTELLKRFISERGKILPRRVTGTSAKNQRKVVNAIKRARVMALLPFVAEDQN
>2_44689
MASTQNIVEEVQKMLDTYDTNKDGEITKAEAVEYFKGKKAFNPERSAIYLFQVYDKDNDGKITIKELAGDIDFDKALKEYKEKQAKSKQQEAEVEEDIEAFILRHNKDDNTDITKDELIQGFKETGAKDPEKSANFILTEMDTNKDGTITVKELRVYYQKVQKLLNPDQ
>3_352472
MKTKSSNNIKKIYYISSILVGIYLCWQIIIQIIFLMDNSIAILEAIGMVVFISVYSLAVAINGWILVGRMKKSSKKAQYEDFYKKMILKSKILLSTIIIVIIVVVVQDIVINFILPQNPQPYVYMIISNFIVGIADSFQMIMVIFVMGELSFKNYFKFKRIEKQKNHIVIGGSSLNSLPVSLPTVKSNESNESNTISINSENNNSKVSTDDTINNVM
>4_91061
MTNPFENDNYTYKVLKNEEGQYSLWPAFLDVPIGWNVVHKEASRNDCLQYVENNWEDLNPKSNQVGKKILVGKR
...
```

:::

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

:::tip
You can speed up database construction by supplying the threads parameter (`-t`).
:::

:::note{title="Expected files in database directory" collapse}

- `kaiju`
  - `kaiju_db_*.fmi`
  - `nodes.dmp`
  - `names.dmp`

:::

For the Kaiju database construction documentation, see the [Kaiju custom database guide](https://github.com/bioinformatics-centre/kaiju#custom-database).

## Kraken2 databases

The Kraken2 database will be used to classify the reads and intermediate contigs in taxonomic groups.

A number of database indexes have already been generated and maintained by [@BenLangmead Lab](https://github.com/BenLangmead); you can browse them in the [AWS Kraken2 index catalogue](https://benlangmead.github.io/aws-indexes/k2). These databases can directly be used to run the workflow with Kraken2 as well as Bracken.

In case the databases above do not contain your desired libraries, you can build a custom Kraken2 database. This requires two components: a taxonomy (consisting of `names.dmp`, `nodes.dmp`, and `*accession2taxid`) files, and the FASTA files you wish to include.

To pull the NCBI taxonomy, you can run the following:

```bash
kraken2-build --download-taxonomy --db <YOUR_DB_NAME>
```

You can then add your FASTA files with the following build command.

```bash
kraken2-build --add-to-library *.fna --db <YOUR_DB_NAME>
```

You can repeat this step multiple times to iteratively add more genomes prior to building.

Once all genomes are added to the library, you can build the database (and optionally clean it up):

```bash
kraken2-build --build --db <YOUR_DB_NAME>
kraken2-build --clean --db <YOUR_DB_NAME>
```

You can then add the `<YOUR_DB_NAME>/` path to your nf-core/taxprofiler database input sheet.

:::tip{title="Expected files in database directory"}

- `kraken2`
  - `opts.k2d`
  - `hash.k2d`
  - `taxo.k2d`

:::

You can follow the Kraken2 [tutorial](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) for a more detailed description.

### Host read removal

nf-core/viralmetagenome uses Kraken2 to remove contaminated reads.

:::info{title="Why kraken2 for host removal?"}
The reason why we use Kraken2 for host removal over regular read mappers is nicely explained in the following papers:

- [Benchmarking of Human Read Removal Strategies for Viral and Microbial Metagenomics](https://www.biorxiv.org/content/10.1101/2025.03.21.644587v1)
- [The human “contaminome”: bacterial, viral, and computational contamination in whole genome sequences from 1000 families](https://www.nature.com/articles/s41598-022-13269-z)
- [Reconstruction of the personal information from human genome reads in gut metagenome sequencing data](https://www.nature.com/articles/s41564-023-01381-3)

:::

The contamination database is likely the largest database. The default databases are made small explicitly to save storage for end users but are not optimal. I would recommend creating a database consisting of the libraries `human, archaea, bacteria` which will be more than 200GB in size. Additionally, it's good practice to include DNA & RNA of the host of origin if known (i.e. mice, ticks, mosquito, ... ). Add it as described above.

Set it with the variable `--host_k2_db`

### Viral Diversity with Kraken2

The metagenomic diversity estimated with Kraken2 is based on the viral refseq database which can cut short if you expect the species within your sample to have a large amount of diversity eg below 80% ANI ([quasi-species](https://link.springer.com/chapter/10.1007/978-3-642-77011-1_1)). To resolve this it's better to create a database that contains a wider species diversity than only one genome per species. Databases that have this wider diversity are [Virosaurus](https://viralzone.expasy.org/8676) or the [RVDB](https://rvdb.dbi.udel.edu/home) which can increase the accuracy of Kraken2. Add it as described above.

Set it with the variable `--kraken2_db`

## Annotation sequences

Identifying the species and the segment of the final genome constructs is done based on a tblastx search (with MMseqs2) to an annotated sequencing dataset. This dataset is by default the [Virosaurus](https://viralzone.expasy.org/8676) as it contains a good representation of the viral genomes and is annotated.

Set the sequences with `--annotation_db`. The annotation data itself comes from a **metadata table** when you pass `--annotation_metadata`, and from the **fasta headers** otherwise.

### Metadata table

`--annotation_metadata` takes a csv or tsv, optionally gzipped. The first column holds the sequence identifiers, every other column becomes a field in the report:

```csv
Accession,species,segment,host
NC_078521,Mykissvirus tructae,3,Oncorhynchus mykiss
KM368312,Alphainfluenzavirus influenzae,3,Sus scrofa
```

Identifiers are matched against the fasta header up to the first whitespace, with or without a version suffix (`NC_078521.1` matches `NC_078521`). Empty cells, and sequences missing from the table, are left unannotated. A column named like one of the search result columns is skipped with a warning, ignoring case - so NCBI Virus' `Length` gives way to the alignment `length` and does not reach the report.

### Fasta headers

Without `--annotation_metadata`, the annotation data has to be in the header itself as `key=value` or `key="value"` pairs, the way Virosaurus ships it:

```text
>754189.6 species="Ungulate tetraparvovirus 3"|segment="nan"|host_common_name="Pig"
>NC_001731; species=Molluscum contagiosum virus; taxid=10279; segment=N/A;
```

- **Quote your values** unless you separate the pairs with `;`. An unquoted value runs to the next `;` or to the end of the line, so `species=Influenza segment=4` becomes a single `species` field.
- Values may not contain a `"` or a `;`, not even inside quotes.
- Keys are single words: `nucleic acid=DNA` is stored as `acid`.
- `key:value` (colon-separated) is not picked up.

### Creating a custom annotation dataset with [BV-BRC](https://www.bv-brc.org/)

[BV-BRC](https://www.bv-brc.org/) carries a lot of metadata and is queried with their [CLI-tool](https://www.bv-brc.org/docs/cli_tutorial/index.html). Here we take every viral reference genome that is not a lab reassortment, metadata and sequences separately:

```bash
# metadata, +/- 5s
p3-all-genomes --eq superkingdom,Viruses --eq reference_genome,Reference --ne host_common_name,'Lab reassortment' --attr genome_id,species,segment,genome_name,genome_length,host_common_name,genbank_accessions,taxon_id > all-virus-anno.txt
# sequences, done separately as it takes much longer to query +/- 1 hour
p3-all-genomes --eq superkingdom,Viruses --eq reference_genome,Reference --ne host_common_name,'Lab reassortment' | p3-get-genome-contigs --attr sequence > all-virus-contigs.txt
```

:::tip
Any attribute can be added to `--attr` and becomes a column in the report. For the complete list see `p3-all-genomes --fields` or their [manual](https://www.bv-brc.org/docs/cli_tutorial/cli_getting_started.html).
:::

`all-virus-anno.txt` is already a tsv whose first column is the `genome_id`, so it can go straight to `--annotation_metadata`. Only the sequences need a conversion step, as `p3-get-genome-contigs` returns a tsv rather than a fasta:

```python
import pandas as pd

contigs = pd.read_csv("all-virus-contigs.txt", sep="\t", dtype=str)
with open("bv-brc-virus.fasta", "w") as out:
    for genome_id, sequence in zip(contigs["genome.genome_id"], contigs["contig.sequence"]):
        out.write(f">{genome_id}\n{sequence}\n")
```

A genome with several contigs is written out once per contig under the same `genome_id`, which is what you want - they share the one metadata row.

The p3 tools prefix their columns with the entity name, so trim that off to keep the report headings readable:

```bash
sed '1s/genome\.//g' all-virus-anno.txt > bv-brc-virus-metadata.tsv
```

```bash
nextflow run nf-core/viralmetagenome \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --annotation_db bv-brc-virus.fasta \
    --annotation_metadata bv-brc-virus-metadata.tsv
```

:::tip{title="Expected files in database directory"}

- `all-virus-contigs.txt` and `all-virus-anno.txt` (the two downloads)
- `bv-brc-virus.fasta` (passed to `--annotation_db`)
- `bv-brc-virus-metadata.tsv` (passed to `--annotation_metadata`)

:::

### Creating a custom annotation dataset with [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)

[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/) is a good alternative to BV-BRC when you want a dataset that is refreshed more often, or that is focused tightly on one taxon of interest. Its FASTA download carries no annotations at all - the defline is just an accession and free text - so the metadata comes from the exported table, which is what `--annotation_metadata` takes.

```text
>NC_078521.1 |Rainbow trout orthomyxovirus-1 strain Rainbow/Idaho/347/1997 putative nucleoprotein (NP) gene, complete cds
```

1. Search "Nucleotide" sequences for your taxon of interest - a species/family name or an NCBI taxid, e.g. [`Orthomyxoviridae, taxid:11308`](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Orthomyxoviridae,taxid:11308).
1. Narrow the result set with the filters (host, sequence length, completeness, RefSeq only, ...) so you end up with a representative, non-redundant set.
1. Enable the metadata columns you want with the column-selector, then `Download` > `Current table` > `CSV` > **Download All Records**. `Accession` comes first and is what the sequences are matched on; `Species`, `Segment` and `Host` describe a consensus genome best.
1. `Download` > `Nucleotide` > `FASTA` > **Download All Records**.
1. Point the pipeline at both files as they are downloaded:

```bash
nextflow run nf-core/viralmetagenome \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --annotation_db sequences.fasta \
    --annotation_metadata sequences.csv
```

:::warning
Both downloads must come from the **same search and the same filters**. Neither file records the query, so a FASTA and a CSV from two different sessions look valid but share no accessions. The pipeline warns how many contigs hit a sequence missing from the metadata table - if that is the whole run, this is why.
:::

:::note{title="What about JenaDB?"}
The best match found for a database going by that name is [VirJenDB](https://www.virjendb.org) (also written VJDB), a FAIR virus (meta)data platform developed at Friedrich Schiller University Jena as part of [NFDI4Microbiota](https://nfdi4microbiota.de/usecases/virjendb) and described in [Nucleic Acids Research](https://academic.oup.com/nar/article/54/D1/D912/8382361). It aggregates and harmonises sequences and roughly 200 metadata fields from 16 upstream sources - including NCBI Virus, BV-BRC, ICTV and ViralZone - and offers FASTA/TSV/CSV/JSON/XML export through its web portal and a documented REST API. It's a comparatively new and still-evolving resource, so check its own documentation for the current export options; the same recipe as above (download the sequences and the metadata table, then pass them as `--annotation_db` and `--annotation_metadata`) applies regardless of which of these sources you pull from.
:::
