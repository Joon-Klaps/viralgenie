#!/usr/bin/env python

"""Provide a command line tool to filter blast results."""

import argparse
import re
import logging
import sys
import pandas as pd
from pathlib import Path

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to extract sequence names from cdhit's cluster files.",
        epilog="Example: python blast_filter.py in.clstr prefix",
    )

    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="cluster file from chdit or vsearch containing cluster information.",
    )

    parser.add_argument(
        "file_out_prefix",
        metavar="FILE_OUT_PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "-e",
        "--escore",
        metavar="escore",
        type=float,
        help="Escore cutoff",
        default=0,
    )

    parser.add_argument(
        "-b",
        "--bitscore",
        metavar="bitscore",
        type=float,
        help="bitscore cutoff",
        default=0,
    )

    parser.add_argument(
        "-a",
        "--percent-alignment",
        metavar="alignment",
        type=float,
        help="percentage of query alignment length cutoff",
        default=0.80,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def filter(df, escore, bitscore, percent_alignment):
    """Filter blast results."""
    if escore != 0:
        df = df[df["evalue"] <= escore]
    if bitscore != 0:
        df = df[df["bitscore"] >= bitscore]
    if percent_alignment != 0:
        df["percent_alignment"] = df["length"] / df["qlen"]
        df = df[df["percent_alignment"] >= percent_alignment]
    return df


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)

    df = pd.read_csv(args.file_in, sep="\t", header=None)
    df.columns = [
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

    df_filter = filter(df, args.escore, args.bitscore, args.percent_alignment)

    df_filter.to_csv(args.file_out_prefix + ".filter.tsv", sep="\t", index=False)
    df_filter["subject"].to_csv(args.file_out_prefix + ".filter.hits.txt", sep="\t", index=False, header=False)


if __name__ == "__main__":
    sys.exit(main())
