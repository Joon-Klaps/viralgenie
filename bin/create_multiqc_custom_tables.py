#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

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
        "--checkv_files",
        metavar="CHECKV FILES",
        nargs="+",
        help="Checkv summary files for each sample",
        type=Path,
    )

    parser.add_argument(
        "--blast_files",
        metavar="BLAST FILES",
        nargs="+",
        help="Blast files for each contig, having the standard outfmt 6",
        type=Path,
    )

    parser.add_argument(
        "--multiqc_dir",
        metavar="MULTIQC DIR",
        help="Multiqc directory where the multiqc files will be used to create the custom tables for multiqc",
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
        if not Path(file).exists():
            logger.error(f"The given input file {file} was not found!")
            sys.exit(2)


def read_header_file(file_path):
    """Read a file and return its content as a list of strings."""
    with open(file_path, "r") as file:
        content = file.read().splitlines()
    return content


def concat_table_files(table_files, **kwargs):
    """Concatenate all the cluster summary files into a single dataframe."""
    df = pd.concat([pd.read_csv(file, sep="\t", **kwargs) for file in table_files])
    return df


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

    # Cluster summaries
    clusters_summary_df = pd.DataFrame()
    if args.clusters_summary:
        header_cluster_summary = f"{args.header_dir}/clusters_summary_mqc.txt"
        # Check if the given files exist
        check_file_exists(args.clusters_summary)
        check_file_exists([header_cluster_summary])

        # Read the header file
        header_cluster_summary = read_header_file(header_cluster_summary)

        # append comments cluster summary to header
        comments_cluster_summary = [
            "# pconfig:",
            "#    id: 'clusters_summary'",
            "#    namespace: 'Clusters summary (" + args.cluster_method + ")'",
            "#    table_title: 'Clusters summary (" + args.cluster_method + ")'",
        ]
        header_cluster_summary.extend(comments_cluster_summary)

        # Concatenate all the cluster summary files into a single dataframe
        clusters_summary_df = concat_table_files(args.clusters_summary)
        write_tsv_file_with_comments(clusters_summary_df, "summary_clusters_mqc.tsv", header_cluster_summary)

    # Sample metadata
    sample_metadata_df = pd.DataFrame()
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

        # Checkv summary
    if args.checkv_files:
        header_checkv = f"{args.header_dir}/checkv_mqc.txt"

        # Check if the given files exist
        check_file_exists(args.checkv_files)
        check_file_exists([header_checkv])

        # Read the header file
        header_checkv = read_header_file(header_checkv)

        # Read the checkv summary file
        checkv_df = concat_table_files(args.checkv_files)

        # Split up the sample names into sample, cluster, step
        checkv_df[["sample", "cluster", "step", "remaining"]] = checkv_df["contig_id"].str.split("_", n=3, expand=True)
        checkv_df.drop(columns=["remaining", "contig_id"], inplace=True)
        checkv_df["step"] = checkv_df["step"].str.split(".").str[0]

        # Reorder the columns
        checkv_df = checkv_df[
            ["sample", "cluster", "step"]
            + [column for column in checkv_df.columns if column not in ["sample", "cluster", "step"]]
        ]

        # Write the dataframe to a file
        write_tsv_file_with_comments(checkv_df, "summary_checkv_mqc.tsv", header_checkv)

    # Blast summary
    if args.blast_files:
        header_blast = f"{args.header_dir}/blast_mqc.txt"

        # Check if the given files exist
        check_file_exists(args.blast_files)
        check_file_exists([header_blast])

        # Read the header file
        header_blast = read_header_file(header_blast)

        # Read the blast summary file
        blast_df = concat_table_files(args.blast_files, header=None)
        blast_df.columns = [
            "query",
            "subject",
            "pident",
            "qlen",
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

        # Filter on best hit per contig and keep only the best hit
        blast_df = blast_df.sort_values("bitscore", ascending=False).drop_duplicates("query")

        # Extract the species name from the subject column
        blast_df["species"] = blast_df["subject"].str.split("|").str[-1]

        #  Split up the sample names into sample, cluster, step
        blast_df[["sample", "cluster", "step", "remaining"]] = blast_df["query"].str.split("_", n=3, expand=True)
        blast_df.drop(columns=["remaining", "query"], inplace=True)
        blast_df["step"] = blast_df["step"].str.split(".").str[0]

        # Reorder the columns
        blast_df = blast_df[
            [
                "sample",
                "cluster",
                "step",
                "species",
                "subject",
                "pident",
                "qlen",
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
        ]

        # Write the dataframe to a file
        write_tsv_file_with_comments(blast_df, "summary_blast_mqc.tsv", header_blast)

    # Multiqc output yml files
    if args.multiqc_dir:
        # Check if the given files exist
        check_file_exists([args.multiqc_dir])

        # Read the multiqc data yml files
        multiqc_data = [file for file in args.multiqc_dir.glob("multiqc_data/*.yml")]

        # Files of interest Sample:
        # TODO: For samples, this is reinventing the wheel, no need to do this, just add to the multiqc general stats table
        files_of_interest = [
            "fastqc",
            "fastp",
            "trimmomatic",
            "bowtie2",
            "bbduk",
            "kraken2",
            "general_stats",
            "sample_metadata",
            "kajiu",
        ]
        # Filter the for the files of interest for samples
        sample_files = [file for file in multiqc_data if any(x in file.stem for x in files_of_interest)]

        # Read in the files
        multiqc_samples_df = pd.DataFrame()
        for file in sample_files:
            with open(file, "r") as stream:
                try:
                    data = yaml.safe_load(stream)
                    df_file = pd.DataFrame()
                    for key, value in data.items():
                        df = pd.DataFrame.from_dict(value)
                        df["sample"] = key
                        df["file"] = file.stem
                        df_file = df_file.concat(df)
                    multiqc_samples_df = multiqc_samples_df.concat(df_file)
                except yaml.YAMLError as exc:
                    print(exc)
        write_tsv_file_with_comments(multiqc_samples_df, "multiqc_output_mqc.tsv", [])

        # # Files of interest contigs:
        # files_of_interest = [
        #     "samtools_stats",
        #     "samtools_idxstats",
        #     "samtools_flagstat",
        # ]

        # Write the dataframe to a file
        # write_tsv_file_with_comments(multiqc_samples_df, "multiqc_output_mqc.tsv", [])

    return 0


if __name__ == "__main__":
    sys.exit(main())
