#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.

"""Provide a command line tool to filter blast results."""

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, cast

import pandas as pd
from Bio import SeqIO
from utils.constant_variables import BLAST_COLUMNS
from utils.file_tools import index_fasta

logger = logging.getLogger()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to filter blast results.",
        epilog="Example: python blast_filter.py -i blast.txt -c contigs.fa -r references.fa -p prefix -b blacklist.txt",
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
        "-k",
        "--blacklist",
        metavar="BLACKLIST FILE",
        type=Path,
        help="File containing subject identifiers to exclude (one per line).",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def filter_hits(
    df: pd.DataFrame,
    escore: float,
    bitscore: float,
    percent_alignment: float,
    blacklist: Optional[set[str]] = None,
) -> pd.DataFrame:
    """Filter blast results."""
    result = df.copy()
    if escore != 0:
        # filter by escore cutoff
        result = cast(pd.DataFrame, result[result["evalue"] <= escore])
    if bitscore != 0:
        # filter by bitscore cutoff
        result = cast(pd.DataFrame, result[result["bitscore"] >= bitscore])
    if percent_alignment != 0:
        # check if alignment percentage meets cutoff
        result["percent_alignment"] = result["length"] / result["qlen"]
        result = cast(pd.DataFrame, result[result["percent_alignment"] >= percent_alignment])
    if blacklist:
        # Compile regex pattern for blacklist entries
        combined_pattern = "|".join([re.escape(entry) for entry in blacklist if entry])
        matches = result["subject"].str.contains(combined_pattern, regex=True)
        result = cast(pd.DataFrame, result[~matches].copy())
        logger.info("Removed %d hits based on blacklist filtering.", matches.sum())
        if result.empty:
            logger.warning("All BLAST hits were removed by blacklist filtering.")
        if combined_pattern and not matches.any():
            logger.debug("Blacklist provided but no valid entries detected after cleaning.")
    return result


def read_blast(blast: Path) -> pd.DataFrame:
    df = pd.read_csv(blast, sep="\t", header=None)
    df.columns = BLAST_COLUMNS
    return df

def load_blacklist(blacklist_path: Path) -> set[str]:
    """Return a set of blacklisted subject identifiers."""
    with open(blacklist_path, "r", encoding="utf-8") as handle:
        entries = {
            line.split()[0]
            for raw_line in handle
            if (line := raw_line.strip()) and not line.startswith("#")
        }
    logger.info("Loaded %d blacklist entries from %s", len(entries), blacklist_path)
    return entries


def write_contigs_and_blast_sequence(
    df: pd.DataFrame,
    contigs: Path,
    references: Path,

    prefix: str,
) -> None:
    """
    Extracts contigs hits from a DataFrame and writes them to a FASTA file.

    Args:
        df (pandas.DataFrame): DataFrame containing the hits information.
        contigs (str): Path to the contigs file.
        references (str): Path to the references file in FASTA format.
        prefix (str): Prefix for the output file.

    Returns:
        None
    """
    # Get unique hit IDs we need to find
    needed_hits = {hit.split(" ")[0] for hit in df["subject"].unique()}

    logger.info("Building sequence index for reference file...")
    ref_index = index_fasta(references)

    # Copy contigs to output file first
    with open(contigs, "r", encoding="utf-8") as contigs_file, open(
        f"{prefix}_withref.fa", "w", encoding="utf-8"
    ) as out_file:
        out_file.write(contigs_file.read())

        # Lookup our sequence
        found_hits: set[str] = set()
        for hit_name in needed_hits:
            if hit_name in ref_index:
                record = ref_index[hit_name]
                # Reassign the id
                record.id = hit_name.replace("|", "-").replace("\\", "-").replace("/", "-")
                SeqIO.write(record, out_file, "fasta")
                found_hits.add(hit_name)
            else:
                logger.debug("Reference sequence %s not found in index", hit_name)

        # Warn if some sequences weren't found
        missing_hits = needed_hits - found_hits
        if missing_hits:
            logger.warning(
                "Could not find the following reference sequences: %s",
                ", ".join(sorted(missing_hits)),
            )

    if hasattr(ref_index, 'close'):
        ref_index.close()


def write_filtered_blast_df(df: pd.DataFrame, prefix: str) -> None:
    # Write filtered hits to file
    df.to_csv(prefix + ".filter.tsv", sep="\t", index=False)

    # Write unique hits to file
    unique_hits = df["subject"].unique()
    unique_series = pd.Series(unique_hits)
    unique_series.to_csv(prefix + ".filter.hits.txt", sep="\t", index=False, header=False)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if args.blast is None and args.contigs and args.contigs.is_file():
        logger.warning("No blast input was provide, just copying input file.")
        with open(args.contigs, "r", encoding="utf-8") as contigs_file:
            contig_content = contigs_file.read()
        with open(f"{args.prefix}_withref.fa", "w", encoding="utf-8") as f:
            f.write(contig_content)
        return 0

    if args.blast is None or not args.blast.is_file():
        logger.error("The given input file %s was not found!", args.blast)
        sys.exit(2)
    if args.references is None or not args.references.is_file():
        logger.error("The given input file %s was not found!", args.references)
        sys.exit(2)
    if args.contigs is None or not args.contigs.is_file():
        logger.error("The given input file %s was not found!", args.contigs)
        sys.exit(2)
    if args.blacklist is not None and not args.blacklist.is_file():
        logger.error("The given blacklist file %s was not found!", args.blacklist)
        sys.exit(2)

    blacklist_hits: Optional[set[str]] = None
    if args.blacklist is not None:
        blacklist_hits = load_blacklist(args.blacklist)

    df = read_blast(args.blast)

    df_filter = filter_hits(
        df,
        args.escore,
        args.bitscore,
        args.percent_alignment,
        blacklist_hits,
    )

    write_contigs_and_blast_sequence(df_filter, args.contigs, args.references, args.prefix)

    write_filtered_blast_df(df_filter, args.prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
