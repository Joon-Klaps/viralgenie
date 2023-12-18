#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract the sequences based on the results of Kaiju, Kraken or both of them combined.",
        epilog="Example: python extract_precluster.py input.tsv sequence.fa prefix",
    )

    parser.add_argument(
        "classifications",
        metavar="CLASSIFICATIONS",
        type=Path,
        help="Classified contigs using Kaiju, Kraken or both of them combined.",
    )

    parser.add_argument(
        "sequences",
        metavar="SEQUENCES",
        type=Path,
        help="Input sequence file used for looking up",
    )

    parser.add_argument(
        "file_out_prefix",
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


def create_groups(file):
    """
    Create groups of sequence names based on classification.

    Args:
        file (str): The path to the input file.

    Returns:
        dict: A dictionary where the keys are taxids and the values are lists of sequence names.
            The 'U' key represents unclassified sequences.
    """
    groups = {}
    with open(file, "r") as f:
        for line in f:
            parts = line.split("\t")
            classified = parts[0]
            seq_name = parts[1]
            taxid = parts[2]
            if classified == "C":
                if taxid not in groups:
                    groups[taxid] = []
                groups[taxid].append(seq_name)
            elif classified == "U":
                groups["U"].append(seq_name)
    return groups


def extract_sequences(groups, sequences, file_out_prefix):
    """
    Extract the sequences from the input file based on the groups.

    Args:
        groups (dict): A dictionary where the keys are taxids and the values are lists of sequence names.
            The 'U' key represents unclassified sequences.
        sequences (str): The path to the input sequence file.
        file_out_prefix (str): The prefix of the output file.
    """
    for taxid, seq_names in groups.items():
        with open(f"{file_out_prefix}_taxid{taxid}.fa", "aw") as f_out:
            for record in SeqIO.parse(sequences, "fasta"):
                if record.id in seq_names:
                    SeqIO.write(record, f_out, "fasta")


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.classifications.is_file():
        logger.error(f"The given input file {args.classifications} was not found!")
        sys.exit(2)
    if not args.sequences.is_file():
        logger.error(f"The given input file {args.sequences} was not found!")
        sys.exit(2)

    groups = create_groups(args.classifications)

    extract_sequences(groups, args.sequences, args.file_out_prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
