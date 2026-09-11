#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.

"""Provide a python file with numerous constant values"""

BLAST_COLUMNS = [
    "query",
    "subject",
    "subject title",
    "pident",
    "qlen",
    "slen",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

READ_DECLARATION = {
    'seqs': {
        'namespace_patterns': 'fastqc',
        'suffix': 'R1,R2'
    },
    'reads after filtering': {
        'namespace_patterns': 'fastp',
        'suffix': 'R1+R2'
    },
    'unique reads': {
        'namespace_patterns': 'humid',
        'suffix': 'R1,R2'
    },
    'reads assigned': {
        'namespace_patterns': 'kaiju',
        'suffix': 'R1+R2'
    }
}


CONSTRAINT_GENERAL_STATS_COLUMNS = [
    "(mash-screen) query-ID",
    "(mash-screen) shared-hashes",
    "(samtools Raw) reads mapped %",
    "(samtools Raw) reads mapped",
    "(samtools Raw) reads unmapped %",
    "(samtools Raw) reads unmapped",
    "(umitools) removed reads",
    "(umitools) deduplicated reads",
    "(picard) read pair duplicates",
    "(samtools Post-dedup) reads mapped %",
    "(samtools Post-dedup) reads mapped",
    "(samtools Post-dedup) reads unmapped %",
    "(samtools Post-dedup) reads unmapped",
    "number_of_SNPs",
    "number_of_indels",
    "qlen",
    "(quast) % N's",
    "CLUSTER: mosdepth.mean_coverage",
    "CLUSTER: mosdepth.min_coverage",
    "CLUSTER: mosdepth.max_coverage",
    "CLUSTER: mosdepth.median_coverage",
    "CLUSTER: mosdepth.1_x_pc",
    "CLUSTER: mosdepth.10_x_pc",
    "CLUSTER: mosdepth.50_x_pc",
    "CLUSTER: mosdepth.100_x_pc",
    "CLUSTER: mosdepth.200_x_pc",
]

COLUMN_MAPPING = {"(blast) qlen": "consensus length", "(annotation) qlen": "consensus length"}

FILES_OF_INTEREST = {
    'samtools-1="samtools Raw"': "",
    "umitools": "multiqc_umitools_dedup",
    "picard": "mutliqc_picard_dups",
    'samtools="samtools Post-dedup"': "multiqc_samtools_stats",
    "ivar_variants": "",
    "bcftools": "multiqc_bcftools_stats",
    "multiqc_general_stats": "",
}

CLUSTER_HEADERS = {
    "# Clusters": {
        "title": "Filtered # clusters",
        "description": "Number of contig clusters used for further refinement",
        "scale": "Blues",
    },
    "# Removed clusters": {
        "title": "Total # clusters",
        "description": "Total number of input contig clusters before filtering",
        "scale": "Blues",
    },
}

CONTIG_TAXONOMY_TOP_N = 5
CONTIG_COMPLETENESS_TOP_N = 10 # Heatmap can fit more

# Species and segment are read back from whatever the annotation database wrote, accomodating for a few
# BV-BRC, virosaurus, and NCBI naming convetions.
ANNOTATION_SPECIES_KEYS = ["species", "speciesname", "virusspecies", "organismname", "organism", "name"]
ANNOTATION_SEGMENT_KEYS = ["segment", "segmentname", "genomesegment", "segmentnumber", "seg"]
ANNOTATION_NULL_VALUES = frozenset({"", "n/a", "na", "nan", "none", "null", "unknown", "undefined", "-", ".", "?"})

# Category names shared by the taxonomy barplot and the completeness heatmap.
TAXON_OTHER = "Other species"
TAXON_UNCLASSIFIED = "Unclassified"
TAXON_UNRECONSTRUCTED = "Not reconstructed"

# Joins a species to its genome segment in the completeness heatmap: "Influenza A virus | seg 4".
SEGMENT_SEPARATOR = "|"

CONTIG_TAXONOMY_PCONFIG = {
    "id": "contig_taxonomy",
    "title": "Contig clusters per sample",
    "ylab": "# clusters",
    "y_decimals": False,
    "use_legend": True,
}

# Sequential (ColorBrewer Greens) instead of the diverging default: pale is an incomplete
# genome, dark is a complete one, which reads more naturally for a 0-100% quality metric.
CONTIG_COMPLETENESS_PCONFIG = {
    "id": "contig_completeness",
    "title": "Consensus genome completeness",
    "xlab": "Species",
    "ylab": "Sample",
    "zlab": "Completeness",
    "min": 0,
    "max": 100,
    "square": False,
    "tt_decimals": 1,
    "cluster_rows": False,
    "cluster_cols": False,
    "xcats_samples": False,
    "colstops": [
        [0, "#f7fcf5"],
        [0.25, "#c7e9c0"],
        [0.5, "#74c476"],
        [0.75, "#31a354"],
        [1, "#006d2c"],
    ],
}


MASH_SCREEN_COLUMNS = ["identity", "shared-hashes", "median-multiplicity", "p-value", "query-ID", "query-comment"]
