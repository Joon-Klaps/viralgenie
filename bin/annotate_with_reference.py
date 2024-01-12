#!/usr/bin/env python

import argparse
import errno
import logging
import sys
from pathlib import Path

import numpy as np
from Bio import Seq, SeqIO, SeqRecord

# global logger
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to annotate the regions with 0 coverage with the reference sequence.",
        epilog="Example: python annotate_with_reference.py --reference reference.fasta --consensus consensus.fasta --mpileup mpileup.txt --prefix prefix",
    )

    parser.add_argument(
        "-r",
        "--reference",
        metavar="REFERENCE FILE",
        type=Path,
        help="Reference sequence file in fasta format.",
    )

    parser.add_argument(
        "-c",
        "--consensus",
        metavar="CONSENSUS FILE",
        type=Path,
        help="Consensus sequence file in fasta format.",
    )

    parser.add_argument(
        "-m",
        "--mpileup",
        metavar="MPILEUP FILE",
        type=Path,
        help="Mpileup file in (default) tsv format.",
    )

    parser.add_argument(
        "-d",
        "--minimum-depth",
        metavar="Minimum depth",
        type=int,
        help="Threshold for minimum depth of coverage (default 0).",
        default=0,
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
        default="WARNING",
    )

    return parser.parse_args(argv)


def annotate_ambiguous(reference_record, consensus_record, regions):
    ambiguous_bases = set("RYKSWMBDHVN")

    ref_seq = np.array(list(str(reference_record.seq)))
    cons_seq = np.array(list(str(consensus_record.seq)))

    mask = (cons_seq == np.array(list(ambiguous_bases))) & np.isin(np.arange(len(cons_seq)), regions)
    cons_seq[mask] = ref_seq[mask]

    return SeqRecord(seq=Seq("".join(cons_seq)), id=consensus_record.id, description="")


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # Read in the reference sequence
    reference = SeqIO.read(args.reference, "fasta")
    # Read in the consensus sequence
    consensus = SeqIO.read(args.consensus, "fasta")

    # Sanity check that the reference and consensus are the same length
    if len(reference) != len(consensus):
        logger.error(f"Reference: {reference} and consensus: {consensus} sequences are not the same length.")
        sys.exit(errno.EINVAL)

    # Read in the mpileup file in a numpy array
    mpileup = np.loadtxt(args.mpileup, dtype=str, delimiter="\t")
    print(mpileup)

    # Extract regions with coverage
    low_coverage = mpileup[mpileup[:, 3].astype(int) <= args.minimum_depth, 1].astype(int)

    print(low_coverage)

    #

    # Create a dictionary of the reference sequence


if __name__ == "__main__":
    sys.exit(main())
