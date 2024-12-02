#!/usr/bin/env python

"""Provide a command line tool to filter generate a custom tsv -csv file from a bam"""

import pysam
import pysamstats
import csv
import argparse
import pandas as pd
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
    parser.add_argument("alignment", type=Path, help="Input BAM file prefix")
    parser.add_argument("reference", type=Path, help="Reference FASTA file")
    parser.add_argument("prefix", type=str, help="Prefix of output file")
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def process_mpileup(filename, reference):
    header = ["Position", "A", "C", "G", "T", "Insertions", "Deletions", "Consensus"]
    if filename.suffix == ".cram":
        alignmentFile = pysam.AlignmentFile(filename, "rc", reference_filename=reference)
    else:
        alignmentFile = pysam.AlignmentFile(filename)  # auto-detects file type

    stats = pysamstats.stat_variation(alignmentFile, fafile=reference)
    logger.debug(list(stats.items())[:10])

    # for record in pysamstats.stat_variation(alignmentFile, fafile=reference):
    #     rec = [record["pos"], record["A"], record["C"], record["G"], record["T"], record["insertions"], record["deletions"]]
    #     if rec[1] + rec[2] + rec[3] + rec[4]:
    #         percentage = float(max(rec[1:5])) / (rec[1] + rec[2] + rec[3] + rec[4])
    #     else:
    #         percentage = 0
    #     ind = rec.index(max(rec[1:5]))
    #     found = ["N", "A", "C", "G", "T"][ind] if ind in [1, 2, 3, 4] else "N"
    #     if max(rec[1:5]) > 19 and percentage >= 0.7:
    #         rec.append(found)
    #         consensus.append(found)
    #     else:
    #         rec.append("N")
    #         consensus.append("N")
    #     results.append(rec)
    # seq = "".join(consensus)

    # with open(f"{csvname}.csv", "w", newline="") as csvfile:
    #     writer = csv.writer(csvfile, delimiter=",")
    #     writer.writerow(header)
    #     writer.writerows(results)

    # with open(f"{fastaname}.fasta", "w") as fastafile:
    #     fastafile.write(f">{filename}\n")
    #     fastafile.write(seq)

    return pd.DataFrame()


def main():
    args = parse_args()
    logger.info("Starting mpileup processing")
    df = process_mpileup(args.filename, args.reference)
    df.to_csv("{prefix}.tsv", sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC)
    logger.info("Mpileup processing completed")


if __name__ == "__main__":
    main()
