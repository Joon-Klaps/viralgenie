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
    parser.add_argument("--k", type=int, help="Pseudocount to add to the total for shannon entropy calculation", default=50)
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def process_mpileup(filename: Path, reference: Path, k: int) -> NDArray:
    """
    Process mpileup data using numpy vectorized operations.
    """
    fasta = pysam.FastaFile(str(reference))
    alignment_file = pysam.AlignmentFile(filename, "rc" if filename.suffix == ".cram" else "rb", reference_filename=str(reference))

    # Convert generator to structured numpy array
    stats = list(pysamstats.stat_variation(alignment_file, fafile=fasta))
    n_rows = len(stats)

    # Create structured array in one go
    data = np.zeros(n_rows, dtype=[
        ("pos", int), ("ref", "U1"), ("A", int), ("C", int),
        ("G", int), ("T", int), ("ins", int), ("del", int),
        ("consensus", "U1"), ("entropy", float), ("weighted_entropy", float)
    ])

    # Fill arrays using vectorized operations
    data["pos"] = np.array([r["pos"] + 1 for r in stats])
    data["ref"] = np.array([r["ref"] for r in stats])
    for base in "ACGT":
        data[base] = np.array([r[base] for r in stats])
    data["ins"] = np.array([r["insertions"] for r in stats])
    data["del"] = np.array([r["deletions"] for r in stats])
    data["consensus"] = "N"

    # Create nucleotide matrix for vectorized operations
    nucleotides = np.stack([data[base] for base in "ACGT"], axis=1)
    total_coverage = np.sum(nucleotides, axis=1)

    # Vectorized consensus calculation
    max_counts = np.max(nucleotides, axis=1)
    mask = np.divide(max_counts, total_coverage, where=total_coverage > 0) >= 0.7
    data["consensus"][mask] = np.array(["A", "C", "G", "T"])[np.argmax(nucleotides[mask], axis=1)]

    # Calculate shannon entropy
    data["entropy"] = shannon_entropy(nucleotides, total_coverage)
    data["weighted_entropy"] = weighted_entropy(data["entropy"], total_coverage, k)

    return data

def shannon_entropy(nucleotides: NDArray, total_coverage: NDArray) -> NDArray:
    """
    Calculate the Shannon entropy of the nucleotide distribution
    """
    # Calculate the frequency of each nucleotide
    frequencies = np.divide(nucleotides, total_coverage[:, np.newaxis], where=total_coverage[:, np.newaxis] > 0)

    # Calculate the Shannon entropy
    with np.errstate(divide='ignore', invalid='ignore'):
        log2_freqs = np.log2(frequencies, where=frequencies > 0)
        entropy = -np.sum(frequencies * log2_freqs, axis=1, where=~np.isnan(log2_freqs))

    # Replace NaN values with 0.0
    entropy = np.nan_to_num(entropy)

    # Replace -0.0 with 0.0
    entropy = np.where(entropy == -0.0, 0.0, np.round(entropy, 3))

    return entropy

def weighted_entropy(entropy: NDArray, total_coverage: NDArray, k: int) -> NDArray:
    """
    Correct the Shannon entropy by multiplying it with N/(N+k)
    """
    correction = total_coverage / (total_coverage + k)
    return np.round(entropy * correction, 3)


def write_csv(matrix: NDArray, prefix: str) -> None:
    """
    Write the matrix to a csv file

    Args:
        matrix: NumPy array containing the mpileup results
        output: Path to the output file
    """
    header = ["Position", "Reference", "A", "C", "G", "T", "Insertions", "Deletions", "Consensus", "Entropy", "Weighted Entropy"]
    with open(f"{prefix}.tsv", "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(header)
        writer.writerows(matrix)


def main():
    args = parse_args()
    logger.info("Starting mpileup processing")
    matrix = process_mpileup(args.alignment, args.reference, args.k)
    write_csv(matrix, args.prefix)
    logger.info("Mpileup processing completed")


if __name__ == "__main__":
    main()
