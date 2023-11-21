#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import logging
import os
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
        "--comment_dir",
        metavar="MULTIQC COMMENT DIR",
        help="Directory with the multiqc header files for table annotation that correspond to the different tables that will be created",
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
        "--quast_files",
        metavar="QUAST FILES",
        nargs="+",
        help="Quast summary files for each sample",
        type=Path,
    )

    parser.add_argument(
        "--save_intermediate",
        metavar="SAVE INTERMEDIATE FILES",
        type=bool,
        nargs="?",
        const=True,
        default=False,
    )

    parser.add_argument(
        "--table_headers",
        metavar="TABLE HEADERS",
        help="Yaml file with the table headers for the different tables that will be created",
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


def check_file_exists(files, throw_error=True):
    """Check if the given files exist."""
    for file in files:
        if not Path(file).exists():
            if throw_error:
                logger.error(f"The given input file {file} was not found!")
                sys.exit(2)
            else:
                logger.warning(f"The given input file {file} was not found!")
        elif not os.stat(file).st_size > 0:
            logger.warning(f"The given input file {file} is empty, it will not be used!")


def read_header_file(file_path):
    """Read a file and return its content as a list of strings."""
    with open(file_path, "r") as file:
        content = file.read().splitlines()
    return content


def concat_table_files(table_files, **kwargs):
    """Concatenate all the cluster summary files into a single dataframe."""
    df = pd.concat([pd.read_csv(file, sep="\t", **kwargs) for file in table_files if os.stat(file).st_size > 0])
    return df


def read_in_quast(table_files):
    """Concatenate all the cluster summary files into a single dataframe."""
    df = pd.DataFrame()
    for file in table_files:
        with open(file, "r") as f:
            d = dict(line.strip().split("\t") for line in f)
        df = pd.concat([df, pd.DataFrame.from_dict(d, orient="index").T])
    return df


def write_tsv_file_with_comments(df, file, comment):
    df_tsv = df.to_csv(sep="\t", index=False)
    with open(file, "w") as f:
        if comment:
            f.write("\n".join(comment))
            f.write("\n")
        f.write(df_tsv)


def read_multiqc_data(directory, files_of_interest):
    # Get all the multiqc data files
    multiqc_data = [file for file in directory.glob("*.txt")]

    # Filter the for the files of interest for contigs
    sample_files = [file for file in multiqc_data if any(x in file.stem for x in files_of_interest)]

    # Read in the files
    multiqc_samples_df = pd.DataFrame()

    # Tsv's
    for file in sample_files:
        with open(file, "r") as table:
            df = pd.read_csv(table, sep="\t")
            df.set_index(df.columns[0], inplace=True)  # Set the first column as index

            if multiqc_samples_df.empty:
                multiqc_samples_df = df
            else:
                multiqc_samples_df = multiqc_samples_df.join(df, how="outer")
    return multiqc_samples_df


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # Cluster summaries
    clusters_summary_df = pd.DataFrame()
    if args.clusters_summary:
        header_cluster_summary = []
        if args.comment_dir:
            header_cluster_summary = f"{args.comment_dir}/clusters_summary_mqc.txt"
            check_file_exists([header_cluster_summary])
            header_cluster_summary = read_header_file(header_cluster_summary)

        # Check if the given files exist
        check_file_exists(args.clusters_summary)

        # Concatenate all the cluster summary files into a single dataframe
        clusters_summary_df = concat_table_files(args.clusters_summary)
        write_tsv_file_with_comments(clusters_summary_df, "summary_clusters_mqc.tsv", header_cluster_summary)

    # Sample metadata
    sample_metadata_df = pd.DataFrame()
    if args.sample_metadata:
        # Check if the given files exist
        check_file_exists([args.sample_metadata])

        header_sample_metadata = []
        if args.comment_dir:
            header_sample_metadata = f"{args.comment_dir}/sample_metadata_mqc.txt"
            check_file_exists([header_sample_metadata])
            header_sample_metadata = read_header_file(header_sample_metadata)

        # Read the sample metadata file
        sample_metadata_df = pd.read_csv(args.sample_metadata, sep="\t")

        # Write the dataframe to a file
        write_tsv_file_with_comments(sample_metadata_df, "sample_metadata_mqc.tsv", header_sample_metadata)

    # Checkv summary
    checkv_df = pd.DataFrame()
    if args.checkv_files:
        header_checkv_list = []
        if args.comment_dir:
            header_checkv = f"{args.comment_dir}/checkv_mqc.txt"
            check_file_exists([header_checkv], throw_error=False)
            try:
                header_checkv_list = read_header_file(header_checkv)
            except:
                header_checkv_list = []

        # Check if the given files exist
        check_file_exists(args.checkv_files)

        # Read the header file
        header_checkv = read_header_file(header_checkv)

        # Read the checkv summary file
        checkv_df = concat_table_files(args.checkv_files)

        # Split up the sample names into sample, cluster, step
        checkv_df = checkv_df.add_prefix("(checkv) ")
        checkv_df[["sample", "cluster", "step", "remaining"]] = checkv_df["(checkv) contig_id"].str.split(
            "_", n=3, expand=True
        )
        checkv_df.drop(columns=["remaining"], inplace=True)
        # Rename for id for joining samples again
        checkv_df["step"] = checkv_df["step"].str.split(".").str[0]
        checkv_df["id"] = checkv_df["sample"] + "_" + checkv_df["cluster"] + "_" + checkv_df["step"]
        checkv_df = checkv_df.set_index("id")

        # Write the dataframe to a file
        if args.save_intermediate:
            write_tsv_file_with_comments(checkv_df, "summary_checkv_mqc.tsv", header_checkv_list)

        checkv_df.drop(columns=["sample", "cluster", "step"], inplace=True)

    # Quast summary
    quast_df = pd.DataFrame()
    if args.quast_files:
        # Check if the given files exist
        check_file_exists(args.quast_files)

        # Read the header file
        header_quast_list = []
        if args.comment_dir:
            header_quast = f"{args.comment_dir}/quast_mqc.txt"
            check_file_exists([header_quast], throw_error=False)
            try:
                header_quast_list = read_header_file(header_quast)
            except:
                header_quast_list = []

        # Read the quast summary files & transpose
        quast_df = read_in_quast(args.quast_files)

        # Split up the sample names into sample, cluster, step
        quast_df = quast_df.add_prefix("(quast) ")
        quast_df[["sample", "cluster", "step"]] = quast_df["(quast) Assembly"].str.split("_", n=2, expand=True)
        quast_df["step"] = quast_df["step"].str.split(".").str[0]
        quast_df["id"] = quast_df["sample"] + "_" + quast_df["cluster"] + "_" + quast_df["step"]
        quast_df = quast_df.set_index("id")

        # Most of the columns are not good for a single contig evaluation
        quast_df = quast_df[["(quast) # N's per 100 kbp"]]

        # Write the dataframe to a file
        if args.save_intermediate:
            write_tsv_file_with_comments(quast_df.reset_index(), "summary_quast_mqc.tsv", header_quast_list)

    # Blast summary
    blast_df = pd.DataFrame()
    if args.blast_files:
        # Check if the given files exist
        check_file_exists(args.blast_files)

        header_blast_list = []
        if args.comment_dir:
            header_blast = f"{args.comment_dir}/blast_mqc.txt"
            check_file_exists([header_blast], throw_error=False)
            try:
                header_blast_list = read_header_file(header_blast)
            except:
                header_blast_list = []

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
        blast_df = blast_df.add_prefix("(blast) ")
        blast_df[["sample", "cluster", "step", "remaining"]] = blast_df["(blast) query"].str.split(
            "_", n=3, expand=True
        )
        blast_df.drop(columns=["remaining", "(blast) query"], inplace=True)
        blast_df["step"] = blast_df["step"].str.split(".").str[0]
        blast_df["id"] = blast_df["sample"] + "_" + blast_df["cluster"] + "_" + blast_df["step"]
        blast_df = blast_df.set_index("id")
        blast_df = blast_df[["(blast) species"] + blast_df.columns.difference(["(blast) species"], sort=False).tolist()]

        # Write the dataframe to a file
        if args.save_intermediate:
            write_tsv_file_with_comments(blast_df, "summary_blast_mqc.tsv", header_blast_list)

        blast_df.drop(columns=["sample", "cluster", "step"], inplace=True)

    # Multiqc output yml files
    if args.multiqc_dir:
        # Check if the given files exist
        check_file_exists([args.multiqc_dir])

        header_clusters_overview = []
        if args.comment_dir:
            header_clusters_overview = f"{args.comment_dir}/contig_overview_mqc.txt"
            check_file_exists([header_clusters_overview], throw_error=True)
            try:
                header_clusters_overview = read_header_file(header_clusters_overview)
            except:
                header_clusters_overview = []

        if args.table_headers:
            check_file_exists([args.table_headers])
            # Read the yaml column annotation file for the different tables
            with open(args.table_headers, "r") as f:
                header_multiqc = yaml.safe_load(f)
            files_of_interest = list(header_multiqc.keys())

            # Flatten, and make dic with old name and new name
            # New name contains the tool and the annotation of the column specified in suppl file
            columns_of_interest = {}
            for bigkey, element in header_multiqc.items():
                if bigkey == "general_stats":
                    bigkey = "multiqc"
                if isinstance(element, list):
                    for item in element:
                        if isinstance(item, dict):
                            for key, value in item.items():
                                columns_of_interest.update({key: f"({bigkey}) {value}"})
                        else:
                            columns_of_interest.update({item: f"({bigkey}) {item.replace('_',' ')}"})
        else:
            # Files of interest contigs:
            files_of_interest = [
                "samtools_stats",
                "umitools",
                "general_stats",
            ]

        # Read the multiqc data yml files
        multiqc_contigs_df = read_multiqc_data(args.multiqc_dir, files_of_interest)

        # Write the complete dataframe to a file
        if args.save_intermediate:
            write_tsv_file_with_comments(multiqc_contigs_df.reset_index(inplace=False), "contigs_intermediate.tsv", [])

        # Filter for columns of interest and rename those we can
        if columns_of_interest:
            columns_available = {
                key: value for key, value in columns_of_interest.items() if key in multiqc_contigs_df.columns
            }
            multiqc_contigs_df = multiqc_contigs_df[columns_available.keys()]
            multiqc_contigs_df.rename(columns=columns_available, inplace=True)

        # Join with the custom contig tables
        multiqc_contigs_df = multiqc_contigs_df.join(checkv_df, how="outer")
        multiqc_contigs_df = multiqc_contigs_df.join(quast_df, how="outer")
        multiqc_contigs_df = multiqc_contigs_df.join(blast_df, how="outer")

        # If we are empty, just quit
        if multiqc_contigs_df.empty:
            logger.warning("No data was found to create the contig overview table!")
            return 0

        # Keep only those rows we can split up in sample, cluster, step
        mqc_contigs_sel = multiqc_contigs_df.reset_index()
        mqc_contigs_sel = mqc_contigs_sel[mqc_contigs_sel["index"].str.contains("_")]
        mqc_contigs_sel[["sample name", "cluster", "step"]] = mqc_contigs_sel["index"].str.split("_", n=3, expand=True)

        # Reorder the columns

        final_columns = (
            ["index", "sample name", "cluster", "step"]
            + [column for column in mqc_contigs_sel.columns if "blast" in column]
            + [column for column in mqc_contigs_sel.columns if "checkv" in column]
            + [column for column in mqc_contigs_sel.columns if "quast" in column]
        )

        remaining_columns = mqc_contigs_sel.columns.difference(final_columns, sort=False).tolist()
        mqc_contigs_sel = mqc_contigs_sel[final_columns + remaining_columns]

        write_tsv_file_with_comments(mqc_contigs_sel, "contigs_overview_mqc.tsv", header_clusters_overview)
    return 0


if __name__ == "__main__":
    sys.exit(main())
