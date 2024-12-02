#!/usr/bin/env python

"""Provide a command line tool to filter generate a custom tsv -csv file from a bam"""

import pysam
import pysamstats
import csv
import numpy as np
from numpy.typing import NDArray
import argparse
from pathlib import Path
import logging

# Initialize logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to create a custom summary of mpileup results from a bam/cram/sam file",
        epilog="Example: python custom_mpileup.py --clusters_summary file1,file2,file3,... ",
    )
    parser = argparse.ArgumentParser(description="Custom mpileup processing script.")
    parser.add_argument("--alignment", type=Path, help="Input BAM file prefix")
    parser.add_argument("--reference", type=Path, help="Reference FASTA file")
    parser.add_argument("--prefix", type=str, help="Name of the output file")
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def process_mpileup(filename: Path, reference: Path) -> NDArray:
    """
    Process mpileup data using numpy vectorized operations.

    Args:
        filename: Path to the alignment file (BAM/CRAM/SAM)
        reference: Path to the reference FASTA file

    Returns:
        NDArray: Array with columns [position, A, C, G, T, insertions, deletions, consensus]
    """
    # Initialize FASTA file properly
    fasta = pysam.FastaFile(str(reference))

    alignment_file = pysam.AlignmentFile(filename, "rc" if filename.suffix == ".cram" else "rb", reference_filename=str(reference))

    # Convert generator to structured numpy array
    data = np.array(
        [
            (r["pos"], r["ref"], r["A"], r["C"], r["G"], r["T"], r["insertions"], r["deletions"], "N")
            for r in pysamstats.stat_variation(alignment_file, fafile=fasta)
        ],
        dtype=[
            ("pos", int),
            ("ref", "U1"),
            ("A", int),
            ("C", int),
            ("G", int),
            ("T", int),
            ("ins", int),
            ("del", int),
            ("consensus", "U1"),
        ],
    )

    # Extract nucleotide counts for consensus calculation
    nucleotides = np.vstack([data[base] for base in "ACGT"]).T
    total_coverage = nucleotides.sum(axis=1)
    max_counts = nucleotides.max(axis=1)

    # Update consensus column where conditions are met
    mask = np.divide(max_counts, total_coverage, where=total_coverage > 0) >= 0.7
    data["consensus"][mask] = np.array(["A", "C", "G", "T"])[nucleotides[mask].argmax(axis=1)]

    return data


def write_csv(matrix: NDArray, prefix: str) -> None:
    """
    Write the matrix to a csv file

    Args:
        matrix: NumPy array containing the mpileup results
        output: Path to the output file
    """
    header = ["Position", "Reference", "A", "C", "G", "T", "Insertions", "Deletions", "Consensus"]
    with open(f"{prefix}.tsv", "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(header)
        writer.writerows(matrix)


def main():
    args = parse_args()
    logger.info("Starting mpileup processing")
    matrix = process_mpileup(args.alignment, args.reference)
    write_csv(matrix, args.prefix)
    logger.info("Mpileup processing completed")


if __name__ == "__main__":
    main()
