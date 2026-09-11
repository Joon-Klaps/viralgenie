#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.

"""Provide a command line tool to filter blast results."""

import argparse
import logging
import sys
import os
from pathlib import Path

import pandas as pd
from utils.constant_variables import MASH_SCREEN_COLUMNS
from utils.file_tools import index_fasta
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


def write_hits(df, references, prefix) -> int:
    """
    Write contigs hits from a DataFrame to a FASTA file using indexed lookups.

    Args:
        df (pandas.DataFrame): DataFrame containing the hits information.
        references (str): Path to the references file in FASTA format.
        prefix (str): Prefix for the output file.

    Returns:
        None
    """
    if df.empty:
        open(f"{prefix}_reference.fa", "a").close()
        open(f"{prefix}.json", "a").close()
        return 0

    needed_hits = set(hit.split(" ")[0] for hit in df["query-ID"].unique())
    found_hits = set()

    logger.info("Building sequence index for reference file...")
    ref_index = index_fasta(references)

    with open(f"{prefix}_reference.fa", "w") as f:
        init_position = f.tell()

        for hit_name in needed_hits:
            if hit_name in ref_index:
                record = ref_index[hit_name]
                # Clean up illegal characters in headers
                record.id = record.id.replace("\\", "-")
                record.description = record.description.replace("\\", "-")
                SeqIO.write(record, f, "fasta")
                found_hits.add(hit_name)

        if f.tell() == init_position:
            logger.error("No reference sequences found in the hits. Exiting...")

        # Warn about missing sequences
        missing_hits = needed_hits - found_hits
        if missing_hits:
            logger.warning(f"Could not find the following reference sequences: {', '.join(missing_hits)}")

    if hasattr(ref_index, 'close'):
        ref_index.close()

    # Writing best hit to JSON for metadata purposes
    df_renamed = df.copy()
    df_renamed["query-ID"] = df_renamed["query-ID"].str.replace("\\", "-")
    df_renamed.to_json(f"{prefix}.json", orient="records", lines=True)

    return 0

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

    if not os.stat(file).st_size > 0:
        return pd.DataFrame()

    try:
        df = pd.read_csv(file, sep="\t", header=None, names=MASH_SCREEN_COLUMNS)
    except pd.errors.EmptyDataError as e:
        logger.warning(f"Empty file: {file}, skipping analysis")
        return pd.DataFrame()

    logger.info("Removing duplicates and sorting by identity and shared-hashes...")
    df["shared-hashes_num"] = df["shared-hashes"].str.split("/").str[0].astype(float)
    df = df.sort_values(by=["identity", "shared-hashes_num"], ascending=False)
    df = df.drop(columns=["shared-hashes_num"])

    return df.iloc[[0]]


def fallback_df(references_file) -> pd.DataFrame:
    """
    Get the ID and description of the first sequence in the references file.

    Args:
        references_file (Path): Path to the references FASTA file

    Returns:
        tuple: (sequence_id, sequence_description) or (None, None) if no sequences found
    """
    record = next(SeqIO.parse(references_file, "fasta"), None)
    if record:
        return pd.DataFrame({
            "identity": [0.0],
            "shared-hashes": ["0/0"],
            "median-multiplicity": [0],
            "p-value": [1.0],
            "query-ID": [record.id],
            "query-comment": [record.description if record.description else ""]
        })
    return pd.DataFrame()


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
        logger.info("No mash screen hits found. Using first reference sequence as fallback.")

        # Get first reference sequence info
        df = fallback_df(args.references)

    # Selecting the best hit and write the hit to a fasta file
    return write_hits(df, args.references, args.prefix)


if __name__ == "__main__":
    sys.exit(main())
