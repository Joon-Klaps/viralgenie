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
    "input_reads",
    "output_reads",
    "number_of_SNPs",
    "CLUSTER: mosdepth.mean_coverage",
    "CLUSTER: mosdepth.min_coverage",
    "CLUSTER: mosdepth.max_coverage",
    "CLUSTER: mosdepth.median_coverage",
    "CLUSTER: mosdepth.1_x_pc",
    "CLUSTER: mosdepth.10_x_pc",
]

FILES_OF_INTEREST = {
    "samtools": "multiqc_samtools_stats",
    "umitools": "multiqc_umitools_dedup",
    "multiqc_general_stats": "",
    "picard": "mutliqc_picard_dups",
    "ivar_variants": "",
    "bcftools": "multiqc_bcftools_stats",
}

CLUSTER_HEADERS = {
    "# Remaining clusters": {
        "title": "Filtered # clusters",
        "description": "Number of contig clusters used for further refinement",
        "scale": "Blues",
    },
    "# Removed clusters": {
        "title": "Total # clusters",
        "description": "Total number of input contig clusters before filtering",
        "scale": "Blues",
    },
    "# Clusters": {
        "title": "# Clusters",
        "description": "Number of contig clusters used for further refinement ",
        "scale": "Blues",
    },
}

CLUSTER_PCONFIG = {
    "id": "summary_clusters_info",
    "title": "Number of contig clusters",
    "ylab": "# clusters",
    "y_decimals": False,
}
