#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger()


def concat_clusters_summary_files(clusters_summary_files):
    """Concatenate all the cluster summary files into a single dataframe."""
    clusters_summary_df = pd.concat([pd.read_csv(file, sep="\t") for file in clusters_summary_files])
    return clusters_summary_df


def write_clusters_summary_file(clusters_summary_df, file, option):
    clusters_summary_tsv = clusters_summary_df.to_csv(sep="\t", index=False)
    with open(f"summary_clusters_mqc.tsv", "w") as f:
        f.write(
            "\n".join(
                [
                    "# id: 'clusters_summary'",
                    "# section_name: 'Clusters summary'",
                    "# format: 'tsv'",
                    "# plot_type: 'table'",
                    "# pconfig:",
                    "#    id: 'clusters_summary'",
                    "#    table_title: 'Clusters summary (" + option + ")'",
                ]
            )
        )
        f.write("\n")
        f.write(clusters_summary_tsv)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to combine individual log & summary files which we will pass down to multiqc subsequently.",
        epilog="Example: python create_multiqc_custom_tables.py --clusters_summary file1,file2,file3,... ",
    )

    parser.add_argument(
        "--clusters_summary",
        metavar="clusters_summary",
        nargs="+",
        type=Path,
        help=" list of cluster summary files from created by the module extract clust.",
    )

    parser.add_argument(
        "--cluster_method",
        metavar="option",
        type=str,
        help=" Algorithm used for clustering of input files",
    )

    parser.add_argument(
        "--prefix",
        metavar="FILE_OUT_PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def check_file_exists(files):
    """Check if the given files exist."""
    for file in files:
        if not file.exists():
            logger.error(f"The given input file {file} was not found!")
            sys.exit(2)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # Cluster summaries
    if args.clusters_summary:
        # Check if the given files exist
        check_file_exists(args.clusters_summary)
        # Concatenate all the cluster summary files into a single dataframe
        clusters_summary_df = concat_clusters_summary_files(args.clusters_summary)
        # Write the dataframe to a file
        write_clusters_summary_file(clusters_summary_df, args.prefix, args.option)
    # Clusters more in depth

    return 0


if __name__ == "__main__":
    sys.exit(main())
