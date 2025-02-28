#!/usr/bin/env python

"""Provide a command line tool to filter blast results."""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from utils.constant_variables import BLAST_COLUMNS

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to filter blast results.",
        epilog="Example: python blast_filter.py in.clstr prefix",
    )

    parser.add_argument(
        "-i",
        "--blast",
        metavar="BLAST FILE",
        type=Path,
        help="Blast result file in specific out format.",
    )

    parser.add_argument(
        "-c",
        "--contigs",
        metavar="CONTIG FILE",
        type=Path,
        help="Contig sequence file that was blasted",
    )

    parser.add_argument(
        "-r",
        "--references",
        metavar="REFERENCE FILE",
        type=Path,
        help="Contig sequence file that was blasted",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        metavar="PREFIX",
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


def read_blast(blast):
    df = pd.read_csv(blast, sep="\t", header=None)
    df.columns = BLAST_COLUMNS
    return df


def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


def extract_contigs_hits(df, contigs, references, prefix):
    """
    Extracts contigs hits from a DataFrame and writes them to a FASTA file.
    Processes reference sequences in chunks to avoid memory issues.

    Args:
        df (pandas.DataFrame): DataFrame containing the hits information.
        contigs (str): Path to the contigs file.
        references (str): Path to the references file in FASTA format.
        prefix (str): Prefix for the output file.

    Returns:
        None
    """
    # Get unique hit IDs we need to find
    needed_hits = set(hit.split(" ")[0] for hit in df["subject"].unique())

    # Copy contigs to output file first
    with open(contigs, "r") as contigs_file, open(f"{prefix}_withref.fa", "w") as out_file:
        out_file.write(contigs_file.read())

        # Process reference sequences in chunks
        found_hits = set()
        for record in SeqIO.parse(references, "fasta"):
            hit_name = record.id
            if hit_name in needed_hits:
                SeqIO.write(record, out_file, "fasta")
                found_hits.add(hit_name)

                # Exit early if we found all needed sequences
                if found_hits == needed_hits:
                    break

        # Warn if some sequences weren't found
        missing_hits = needed_hits - found_hits
        if missing_hits:
            logger.warning(f"Could not find the following reference sequences: {', '.join(missing_hits)}")


def write_hits(df, contigs, references, prefix):
    # Extract contigs & blast hits and write to fasta file
    extract_contigs_hits(df, contigs, references, prefix)

    # Write filtered hits to file
    df.to_csv(prefix + ".filter.tsv", sep="\t", index=False)

    # Write unique hits to file
    unique_hits = df["subject"].unique()
    unique_series = pd.Series(unique_hits)
    unique_series.to_csv(prefix + ".filter.hits.txt", sep="\t", index=False, header=False)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if args.blast is None and args.contigs.is_file():
        logger.warning(f"No blast input was provide, just copying input file.")
        with open(args.contigs, "r") as contigs_file:
            contig_content = contigs_file.read()
        with open(f"{args.prefix}_withref.fa", "w") as f:
            f.write(contig_content)
        return 0

    if not args.blast.is_file():
        logger.error(f"The given input file {args.blast} was not found!")
        sys.exit(2)
    if not args.references.is_file():
        logger.error(f"The given input file {args.references} was not found!")
        sys.exit(2)
    if not args.contigs.is_file():
        logger.error(f"The given input file {args.contigs} was not found!")
        sys.exit(2)

    df = read_blast(args.blast)

    df_filter = filter(df, args.escore, args.bitscore, args.percent_alignment)

    write_hits(df_filter, args.contigs, args.references, args.prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
