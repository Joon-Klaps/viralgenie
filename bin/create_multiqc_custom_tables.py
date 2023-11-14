#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to combine individual log & summary files which we will pass down to multiqc subsequently.",
        epilog="Example: python create_multiqc_custom_tables.py --clusters_summary file1,file2,file3,... ",
    )

    parser.add_argument(
        "--clusters_summary",
        metavar="CLUSTER SUMMARY FILES",
        nargs="+",
        type=Path,
        help=" list of cluster summary files from created by the module extract clust.",
    )

    parser.add_argument(
        "--cluster_method",
        metavar="METHOD",
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
        "--sample_metadata",
        metavar="SAMPLE METADATA",
        help="Sample metadata file",
        type=Path,
    )

    parser.add_argument(
        "--header_dir",
        metavar="MULTIQC HEADER DIR",
        help="Directory with the multiqc header files that correspond to the different tables that will be created",
        type=Path,
    )

    parser.add_argument(
        "--sample_metadata",
        metavar="SAMPLE METADATA",
        help="Sample metadata file",
        type=Path,
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


def read_header_file(file_path):
    """Read a file and return its content as a list of strings."""
    with open(file_path, "r") as file:
        content = file.readlines()
    return content


def concat_clusters_summary_files(clusters_summary_files):
    """Concatenate all the cluster summary files into a single dataframe."""
    clusters_summary_df = pd.concat([pd.read_csv(file, sep="\t") for file in clusters_summary_files])
    return clusters_summary_df


def write_tsv_file_with_comments(df, file, comment):
    df_tsv = df.to_csv(sep="\t", index=False)
    with open(file, "w") as f:
        f.write("\n".join(comment))
        f.write("\n")
        f.write(df_tsv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    headers = os.listdir(args.header_dir)

    # Cluster summaries
    if args.clusters_summary:
        header_cluster_summary = f"{args.header_dir}/clusters_summary_mqc.txt"
        # Check if the given files exist
        check_file_exists(args.clusters_summary)
        check_file_exists([header_cluster_summary])

        # Read the header file
        header_cluster_summary = read_header_file(header_cluster_summary)

        # append comments cluster summary to header
        comments_cluster_summary = [
            "#    id: 'clusters_summary'",
            "#    namespace: 'Clusters summary (" + args.cluster_method + ")'",
            "#    table_title: 'Clusters summary (" + args.cluster_method + ")'",
        ]
        header_cluster_summary.extend(comments_cluster_summary)

        # Concatenate all the cluster summary files into a single dataframe
        clusters_summary_df = concat_clusters_summary_files(args.clusters_summary)
        write_tsv_file_with_comments(clusters_summary_df, "summary_clusters_mqc.tsv", comments_cluster_summary)

    # Sample metadata
    if args.sample_metadata:
        header_sample_metadata = f"{args.header_dir}/sample_metadata_mqc.txt"

        # Check if the given files exist
        check_file_exists([args.sample_metadata])
        check_file_exists([header_sample_metadata])

        # Read the header file
        header_sample_metadata = read_header_file(header_sample_metadata)

        # Read the sample metadata file
        sample_metadata_df = pd.read_csv(args.sample_metadata, sep="\t")

        # Write the dataframe to a file
        write_tsv_file_with_comments(sample_metadata_df, "sample_metadata_mqc.tsv", header_sample_metadata)

    # Clusters more in depth
    return 0


if __name__ == "__main__":
    sys.exit(main())
