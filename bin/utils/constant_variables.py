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

CONSTRAIN_GENERAL_STATS_COLUMNS = [
    "(samtools) reads mapped",
    "(samtools) reads_unmapped",
    "number_of_SNPs",
    "number_of_indels" "CLUSTER: mosdepth.mean_coverage",
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
]

COLUMN_MAPPING = {"(blast) qlen": "consensus length", "(annotation) qlen": "consensus length"}

FILES_OF_INTEREST = {
    "samtools": "multiqc_samtools_stats",
    "umitools": "multiqc_umitools_dedup",
    "multiqc_general_stats": "",
    "picard": "mutliqc_picard_dups",
    "ivar_variants": "",
    "bcftools": "multiqc_bcftools_stats",
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
