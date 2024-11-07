#!/usr/bin/env python

"""Provide a command line tool to filter blast results."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger()


def parse_args(argv=None) -> argparse.Namespace:
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to filter blast results.",
        epilog="Example: python blast_filter.py in.clstr prefix",
    )

    parser.add_argument(
        "-i",
        "--mash",
        metavar="MASH FILE",
        type=Path,
        help="Mash screen result file.",
    )

    parser.add_argument(
        "-r",
        "--references",
        metavar="REFERENCE FILE",
        type=Path,
        help="Contig sequence file was screened against",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        metavar="PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def to_dict_remove_dups(sequences) -> dict:
    return {record.id: record for record in sequences}


def write_hits(df, references, prefix) -> None:
    """
    write contigs hits from a DataFrame and writes them to a FASTA file.

    Args:
        df (pandas.DataFrame): DataFrame containing the hits information.
        contigs (str): Path to the contigs file.
        references (str): Path to the references file in FASTA format.
        prefix (str): Prefix for the output file.

    Returns:
        None
    """
    try:
        ref_records = SeqIO.to_dict(SeqIO.parse(references, "fasta"))
    except ValueError as e:
        logger.warning(
            "Indexing the reference pool file causes an error: %s \n Make sure all fasta headers are unique and it is in fasta format! \n AUTOFIX: Taking last occurence of duplicates to continue analysis",
            e,
        )
        ref_records = to_dict_remove_dups(SeqIO.parse(references, "fasta"))
    with open(f"{prefix}_reference.fa", "w") as f:
        init_position = f.tell()
        for hit in df["query-ID"].unique():
            hit_name = hit.split(" ")[0]
            if hit_name in ref_records:
                # Sometimes reads can have illegal characters in the header
                ref_records[hit_name].id = ref_records[hit_name].id.replace("\\", "_")
                ref_records[hit_name].description = ref_records[hit_name].description.replace("\\", "_")
                SeqIO.write(ref_records[hit_name], f, "fasta")
        if f.tell() == init_position:
            logger.error("No reference sequences found in the hits. Exiting...")


def read_mash_screen(file) -> pd.DataFrame:
    """
    Read in the file and return a pandas DataFrame
    File format:
    [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
    0.996341	3786/4000	121	0	USANYPRL230901_81A112023
    0.997144	3832/4000	124	0	EnglandQEUH3267E6482022
    0.997039	3826/4000	121	0	OP958840
    0.997022	3825/4000	122	0	OP971202
    """

    logger.info("Reading in the mash screen file...")
    try:
        df = pd.read_csv(file, sep="\t", header=None)
    except pd.errors.EmptyDataError as e:
        logger.warning(f"Empty file: {file}, skipping analysis")
        return pd.DataFrame()

    df.columns = ["identity", "shared-hashes", "median-multiplicity", "p-value", "query-ID", "query-comment"]

    logger.info("Removing duplicates and sorting by identity and shared-hashes...")
    df["shared-hashes_num"] = df["shared-hashes"].str.split("/").str[0].astype(float)
    df = df.sort_values(by=["identity", "shared-hashes_num"], ascending=False)
    df = df.drop(columns=["shared-hashes_num"])

    return df.iloc[[0]]


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.mash.is_file():
        logger.error(f"The given input file {args.mash} was not found!")
        sys.exit(2)
    if not args.references.is_file():
        logger.error(f"The given input file {args.references} was not found!")
        sys.exit(2)

    # reading in the mash results
    df = read_mash_screen(args.mash)
    if df.empty:
        # Create empty files so nextflow doesn't crash
        open(f"{args.prefix}_reference.fa", "a").close()
        open(f"{args.prefix}.json", "a").close()
        return 0

    # Selecting the best hit and write the hit to a fasta file
    write_hits(df, args.references, args.prefix)

    # Writing the best hit to a json file for metadata purposes
    df.to_json(f"{args.prefix}.json", orient="records", lines=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
