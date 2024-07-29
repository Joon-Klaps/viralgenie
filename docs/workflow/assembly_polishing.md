
# Assembly & polishing

Viralgenie offers an elaborate workflow for the assembly and polishing of viral genomes:

- [Assembly](#de-novo-assembly): combining the results of multiple assemblers.
- [Reference Matching](#reference-matching): comparing the newly assembled contigs to a reference sequence pool.
- [Clustering](#clustering): clustering the contigs based on taxonomy and similarity.
- [Scaffolding](#scaffolding): scaffolding the contigs to the centroid of each bin.
- [Annotation with Reference](#annotation-with-reference): annotating regions with 0-depth coverage with the reference sequence.

![assembly_polishing](../images/assembly_polishing.png)

> The overal workflow of creating reference assisted assemblies can be skipped with the argument `--skip_assembly`. See the [parameters assembly section](../parameters.md#assembly) for all relevant arguments to control the assembly steps.

> The overall refinement of contigs can be skipped with the argument `--skip_polishing`. See the [parameters polishing section](../parameters.md#polishing) for all relevant arguments to control the polishing steps.

The consensus genome of all clusters are then send to the [variant analysis & iterative refinement](variant_and_refinement.md) step.

## De-novo Assembly
Three assemblers are used, [SPAdes](http://cab.spbu.ru/software/spades/), [Megahit](https://github.com/voutcn/megahit), and [Trinity](https://github.com/trinityrnaseq/trinityrnaseq). The resulting contigs of all specified assemblers, are combined and processed further together.
> Modify the spades mode with `--spades_mode [default: rnaviral]` and supply specific params with `--spades_yml` or a hmm model with `--spades_hmm`.

> Specify the assemblers to use with the `--assemblers` parameter where the assemblers are separated with a ','. The default is `spades,megahit,trinity`.

Contigs can be extended using [SSPACE Basic](https://github.com/nsoranzo/sspace_basic) with the `--skip_sspace_basic false` parameter. SSPACE is a tool for scaffolding contigs using paired-end reads. It is modified from SSAKE assembler and has the feature of extending contigs using reads that are unmappable in the contig assembly step.

Low complexity contigs can be filtered out using prinseq++ with the `--skip_contig_prinseq false` parameter. Complexity filtering is primarily a run-time optimisation step. Low-complexity sequences are defined as having commonly found stretches of nucleotides with limited information content (e.g. the dinucleotide repeat CACACACACA). Such sequences can produce a large number of high-scoring but biologically insignificant results in database searches. Removing these reads therefore saves computational time and resources.

## Reference Matching
The newly assembled contigs are compared to a reference sequence pool (--reference_pool) using a [BLASTn search](https://www.ncbi.nlm.nih.gov/books/NBK153387/). This process not only helps annotate the contigs but also assists in linking together sets of contigs that are distant within a single genome. Essentially, it aids in identifying contigs belonging to the same genomic segment and choosing the right reference for scaffolding purposes.

The top 5 hits for each contig are combined with the denovo contigs and send to the clustering step.

> The reference pool can be specified with the `--reference_pool` parameter. The default is the latest clustered [Reference Viral DataBase (RVDB)](https://rvdb.dbi.udel.edu/).

## Clustering

The clustering workflow of contigs consists out of 2 steps, the [pre-clustering](#pre-clustering) and [actual clustering](#actual-clustering). Here contigs are first separated based on identified taxonomy-id ([Kraken2](https://ccb.jhu.edu/software/kraken2/), [Kaiju](https://kaiju.binf.ku.dk/)) and are subsequently clustered further to identify genome segments.

### Pre-clustering

The contigs along with their references have their taxonomy assigned using [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Kaiju](https://kaiju.binf.ku.dk/).

> The default databases are the same ones used for read classification:
> - Kraken2: viral refseq database, `--kraken2_db`
> - Kaiju: clustered [RVDB](https://rvdb.dbi.udel.edu/), `--kaiju_db`

As Kajiu and Kraken2 can have different taxonomic assignments, an additional step is performed to resolve potential inconsistencies in taxonomy and to identify the taxonomy of the contigs. This is done with a custom script that is based on `KrakenTools extract_kraken_reads.py` and `kaiju-Merge-Outputs`.

```mermaid
graph LR;
    A[Contigs] --> B["`**Kraken2**`"];
    A --> C["`**Kaiju**`"];
    B --> D[Taxon merge resolving];
    C --> D;
    D --> E["Taxon filtering"];
    E --> F["Taxon simplification"];
```

!!! Tip annotate "Having very complex metagenomes"
    The pre-clustering step can be used to simplify the taxonomy of the contigs, let [NCBI's taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) help you identify taxon-id's for simplification. The simplification can be done in several ways:

    - Make sure your contamination database is up to date and removes the relevant taxa.
    - Exclude unclassified contigs with `--keep_unclassified false` parameter.
    - Simplify the taxonomy of the contigs to a higher rank using `--precluster_simplify_taxa` parameter (1).
    - Specify the taxa to include or exclude with `--precluster_include_children`(2), `--precluster_include_parents`(3), `--precluster_exclude_children`, `--precluster_exclude_parents`, `--precluster_exclude_taxa` parameters.
    !!! warning
        Providing lists to nextflow is done by encapsulating values with `"` and separating them with a space. For example: `--precluster_exclude_taxa "taxon1 taxon2 taxon3"`.

1. Options here are 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom' or 'superkingdom'.

2. `--precluster_include_childeren`__"genus1"__ :

    ```mermaid
    graph TD;
        A[family] -.- B["genus1 (included)"];
        A -.- C[genus2];
        B --- D[species1];
        B --- E[species2];
        C -.- F[species3];
    ```

3. `--precluster_include_parents` __"species3"__ :

    ```mermaid
    graph TD;
        A["family (included)"] -.- B["genus1"]
        A --- C[genus2]
        B -.- D[species1]
        B -.- E[species2]
        C --- F[species3]
    ```

> The pre-clustering step will be run by default but can be skipped with the argument `--skip_preclustering`. Specify which classifier to use with `--precluster_classifiers` parameter. The default is `kaiju,kraken2`. Contig taxon filtering is still enabled despite not having to solve for inconsistencies if only Kaiju or Kraken2 is ran.

### Actual clustering

The clustering is performed with one of the following tools:

- [`cdhitest`](https://sites.google.com/view/cd-hit)
- [`vsearch`](https://github.com/torognes/vsearch/wiki/Clustering)
- [`mmseqs-linclust`](https://github.com/soedinglab/MMseqs2/wiki#linear-time-clustering-using-mmseqs-linclust)
- [`mmseqs-cluster`](https://github.com/soedinglab/MMseqs2/wiki#cascaded-clustering)
- [`vRhyme`](https://github.com/AnantharamanLab/vRhyme)
- [`mash`](https://github.com/marbl/Mash)


These methods all come with their own advantages and disadvantages. For example, cdhitest is very fast but cannot be used for large viruses >10Mb and similarity threshold cannot go below 80% which is not preferable for highly diverse RNA viruses. Vsearch is slower but accurate. Mmseqs-linclust is the fastest but tends to create a large amount of bins. Mmseqs-cluster is slower but can handle larger datasets and is more accurate. vRhyme is a new method that is still under development but has shown promising results but can sometimes not output any bins when segments are small. Mash is a very fast comparison method is linked with a custom script that identifies communities within a network.

!!! Tip
    When pre-clustering is performed, it is recommended to set a lower identity_threshold (60-70% ANI) as the new goal becomes to separate genome segments within the same bin.

> The clustering method can be specified with the `--clustering_method` parameter. The default is `mash`.

> The network clustering method for `mash` can be specified with the `--network_clustering` parameter. The default is `connected_components`, alternative is [`leiden`](https://www.nature.com/articles/s41598-019-41695-z).

> The similarity threshold can be specified with the `--similarity_threshold` parameter. The default is `0.6`. However, for cdhit the default is `0.8` which is its minimum value.

## Scaffolding

After classifying all contigs and their top BLAST hits into distinct clusters or bins, the contigs are then scaffolded to the centroid of each bin. Any external references that are not centroids of the cluster are subsequently removed to prevent further bias. All members of the cluster are consequently mapped towards their centroid with [Minimap2](https://github.com/lh3/minimap2) and conensus is called using [iVar-consensus](https://andersen-lab.github.io/ivar/html/manualpage.html).


## Annotation with Reference

Regions with 0-depth coverage are annotated with the reference sequence. This is done with a [custom script](https://github.com/Joon-Klaps/viralgenie/blob/dev/bin/lowcov_to_reference.py) that uses the coverage of the denovo contigs towards the reference sequence to identify regions with 0-depth coverage. The reference sequence is then annotated to these regions.

> This step can be skipped using `--skip_hybrid_consensus` parameter.
