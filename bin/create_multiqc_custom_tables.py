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


def check_file_exists(file, throw_error=True):
    """Check if the given files exist."""
    if not Path(file).exists():
        if throw_error:
            logger.error(f"The given input file {file} was not found!")
            sys.exit(2)
        else:
            logger.warning(f"The given input file {file} was not found!")
            return False
    elif not os.stat(file).st_size > 0:
        logger.warning(f"The given input file {file} is empty, it will not be used!")
        return False
    return True


def read_header_file(file_path):
    """Read a file and return its content as a list of strings."""
    with open(file_path, "r") as file:
        content = file.read().splitlines()
    return content


def concat_table_files(table_files, **kwargs):
    """Concatenate all the cluster summary files into a single dataframe."""
    df = pd.concat([pd.read_csv(file, sep="\t", **kwargs) for file in table_files if check_file_exists(file)])
    return df


def read_in_quast(table_files):
    """Concatenate all the cluster summary files into a single dataframe."""
    df = pd.DataFrame()
    for file in table_files:
        d = {}
        if check_file_exists(file, throw_error=False):
            with open(file, "r") as f:
                d = dict(line.strip().split("\t") for line in f)
        df = pd.concat([df, pd.DataFrame.from_dict(d, orient="index").T])
    return df


def write_dataframe(df, file, comment):
    df_tsv = df.to_csv(sep="\t", index=False)
    with open(file, "w") as f:
        if comment:
            f.write("\n".join(comment))
            f.write("\n")
        f.write(df_tsv)


def read_multiqc_data(directory, files_of_interest):
    # Get all the multiqc data files
    multiqc_data = [file for file in directory.glob("multiqc_*.txt")]

    # Filter the for the files of interest for contigs
    sample_files = [file for file in multiqc_data if any(x in file.stem for x in files_of_interest)]

    # Read in the files
    multiqc_samples_df = pd.DataFrame()

    # Tsv's
    for file in sample_files:
        with open(file, "r") as table:
            df = pd.read_csv(table, sep="\t")
            df[df.columns[0]] = df[df.columns[0]].str.split(".").str[0]
            df.set_index(df.columns[0], inplace=True)  # Set the first column as index

            if multiqc_samples_df.empty:
                multiqc_samples_df = df
            else:
                multiqc_samples_df = multiqc_samples_df.join(df, how="outer")
    return multiqc_samples_df


def read_failed_contig_data(directory, files_of_interest):
    multiqc_data = [file for file in directory.glob("multiqc_*.txt")]

    # Filter the for the files of interest for contigs
    sample_files = [file for file in multiqc_data if any(x in file.stem for x in files_of_interest)]

    # Read in the files
    multiqc_samples_df = pd.DataFrame()

    # Tsv's
    for file in sample_files:
        with open(file, "r") as table:
            df = pd.read_csv(table, sep="\t")
            df["Cluster"] = df["Cluster"].str.split(".0").str[0]
            df["Iteration"] = df["Iteration"].str.split(".0").str[0]
            df["id"] = df["Sample"] + "_" + df["Cluster"] + "_" + df["Iteration"]
            df.set_index("id", inplace=True)  # Set the first column as index
            df = df[["id"]]

            if multiqc_samples_df.empty:
                multiqc_samples_df = df
            else:
                multiqc_samples_df = multiqc_samples_df.join(df, how="outer")

    dic = {"Contig qc fail": "True"}
    multiqc_samples_df.assign(**dic)
    return


def get_header(comment_dir, header_file_name):
    header = []
    if comment_dir:
        header_file_path = f"{comment_dir}/{header_file_name}"
        if check_file_exists(header_file_path):
            header = read_header_file(header_file_path)
    return header


# Function to split column dynamically
def dynamic_split(row, delimiter="_"):
    parts = row.split(delimiter)
    return parts


def handle_tables(table_files, header_name=False, output=False, **kwargs):
    result_df = pd.DataFrame()
    if table_files:
        result_df = concat_table_files(table_files, **kwargs)
    if output:
        write_dataframe(result_df, output, header_name)
    return result_df


def handle_dataframe(df, prefix, column_to_split, header=False, output=False):
    result_df = pd.DataFrame()
    if not df.empty:
        df = df.add_prefix(f"({prefix}) ")
        # Apply the dynamic split function to each row in the column
        df[[f"part_{i}" for i in range(1, 5)]] = (
            df[f"({prefix}) {column_to_split}"].apply(dynamic_split, delimiter="_").apply(pd.Series)
        )
        df.rename(
            columns={"part_1": "sample", "part_2": "cluster", "part_3": "step", "part_4": "remaining"}, inplace=True
        )

        df.drop(columns=["remaining"], inplace=True)

        df["step"] = df["step"].str.split(".").str[0]
        df["id"] = df["sample"] + "_" + df["cluster"] + "_" + df["step"]
        df = df.set_index("id")
        result_df = df
    if output:
        write_dataframe(result_df, output, [])
    result_df.drop(columns=["sample", "cluster", "step"], inplace=True)
    return result_df


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # Cluster summaries
    cluster_header = get_header(args.comment_dir, "clusters_summary_mqc.txt")
    clusters_summary_df = handle_table([args.clusters_summary], cluster_header, "summary_clusters_mqc.tsv")

    # Sample metadata
    sample_header = get_header(args.comment_dir, "sample_metadata_mqc.txt")
    sample_metadata_df = handle_table([args.sample_metadata], sample_header, "sample_metadata_mqc.tsv")

    # Checkv summary
    checkv_df = handle_table([args.checkv_files])
    checkv_header = []
    if args.save_intermediate:
        checkv_header = get_header(args.comment_dir, "checkv_mqc.txt")
    if not checkv_df.empty:
        checkv_df = handle_dataframe(checkv_df, "checkv", "contig_id", checkv_header, "summary_checkv_mqc.tsv")

    # Quast summary
    quast_df = read_in_quast([args.quast_files])
    quast_header = []
    if args.save_intermediate:
        quast_header = get_header(args.comment_dir, "quast_mqc.txt")
    if not quast_df.empty:
        quast_df = handle_dataframe(quast_df, "quast", "Assembly", quast_header, "summary_quast_mqc.tsv")

        # Most of the columns are not good for a single contig evaluation
        quast_df = quast_df[["(quast) # N's per 100 kbp"]]

    # Blast summary
    blast_df = handle_table([args.blast_files], header=None)
    blast_header = []
    if args.save_intermediate:
        blast_header = get_header(args.comment_dir, "blast_mqc.txt")
    if not blast_df.empty:
        # Read the blast summary file
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
        blast_df = blast_df[["species"] + blast_df.columns.difference(["species"], sort=False).tolist()]

        blast_df = handle_dataframe(blast_df, "blast", "query", blast_header, "summary_blast_mqc.tsv")

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
                "picard_dups",
                "ivar_variants",
                "bcftools_stats",
            ]

        # Read the multiqc data yml files
        multiqc_contigs_df = read_multiqc_data(args.multiqc_dir, files_of_interest)

        # If we are empty, just quit
        if multiqc_contigs_df.empty:
            logger.warning("No data was found to create the contig overview table!")
            return 0

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

        # adding a tag saying that contig faild qc check
        failed_contigs = ["failed_mapped", "failed_contig_quality"]
        failed_contigs_df = read_failed_contig_data(args.multiqc_dir, failed_contigs)

        # If we are empty, just quit
        if multiqc_contigs_df.empty:
            logger.warning("No data was found to create the contig overview table!")
            return 0

        # Keep only those rows we can split up in sample, cluster, step
        mqc_contigs_sel = multiqc_contigs_df.reset_index().rename(columns={multiqc_contigs_df.index.name: "index"})
        mqc_contigs_sel = mqc_contigs_sel[mqc_contigs_sel["index"].str.contains("_", na=False)]
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
