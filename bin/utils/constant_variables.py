#!/usr/bin/env python

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

}


CONSTRAINT_GENERAL_STATS_COLUMNS = [
    "(samtools Post-dedup) reads mapped %",
    "(samtools Post-dedup) reads mapped",
    "(samtools Post-dedup) reads unmapped %",
    "(samtools Post-dedup) reads unmapped",
    "(samtools Raw) reads mapped %",
    "(samtools Raw) reads mapped",
    "(samtools Raw) reads unmapped %",
    "(samtools Raw) reads unmapped",
    "number_of_SNPs",
    "number_of_indels",
    "CLUSTER: mosdepth.mean_coverage",
    "CLUSTER: mosdepth.min_coverage",
    "CLUSTER: mosdepth.max_coverage",
    "CLUSTER: mosdepth.median_coverage",
    "CLUSTER: mosdepth.1_x_pc",
    "CLUSTER: mosdepth.10_x_pc",
    "CLUSTER: mosdepth.50_x_pc",
    "CLUSTER: mosdepth.100_x_pc",
    "CLUSTER: mosdepth.200_x_pc",
    "qlen",
    "(quast) % N's",
    "(mash-screen) query-ID",
    "(mash-screen) shared-hashes",
    "(failed_mapped) mapped reads",
    "(umitools) removed reads",
    "(umitools) deduplicated reads"
]

COLUMN_MAPPING = {"(blast) qlen": "consensus length", "(annotation) qlen": "consensus length"}

FILES_OF_INTEREST = {
    "failed_mapped": "",
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

CLUSTER_PCONFIG = {
    "id": "summary_clusters_info",
    "title": "Number of contig clusters",
    "ylab": "# clusters",
    "y_decimals": False,
}


MASH_SCREEN_COLUMNS = ["identity", "shared-hashes", "median-multiplicity", "p-value", "query-ID", "query-comment"]
