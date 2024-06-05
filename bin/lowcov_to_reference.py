#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# global logger
logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to annotate the regions with 0 coverage with the reference sequence.",
        epilog="Example: python lowcov_to_reference.py --reference reference.fasta --consensus consensus.fasta --mpileup mpileup.txt --prefix prefix",
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
        help="Mpileup file in (default) tsv format, typically from iVar consensus.",
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
        default="INFO",
    )

    return parser.parse_args(argv)


def annotate_ambiguous(reference, consensus, regions, args):
    """
    Annotates ambiguous bases in the consensus sequence with the corresponding bases from the reference sequence.

    Args:
        reference (SeqRecord): The reference sequence record.
        consensus (SeqRecord): The consensus sequence record.
        regions (list): A list of indices representing the regions where annotation should be performed.

    Returns:
        SeqRecord: The annotated consensus sequence record.
    """
    ID = f"{args.prefix}"
    DESCRIPTION = f"Hybrid construct of the {args.reference} and {args.consensus} sequences where regions with depth lower then {args.minimum_depth} have been replaced"

    seq = alignment_replacement(reference, consensus, regions)

    return SeqRecord(seq=Seq(seq), id=ID, description=DESCRIPTION)

def alignment_replacement(reference_record, consensus_record, regions):
    """
    Replaces the aligned regions in the consensus sequence with the corresponding regions from the reference sequence.

    Args:
        reference_record (Bio.SeqRecord.SeqRecord): The reference sequence record.
        consensus_record (Bio.SeqRecord.SeqRecord): The consensus sequence record.
        regions (list): List of regions to be replaced.

    Returns:
        str: The updated consensus sequence.
    """
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1.0
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -7.0
    aligner.extend_gap_score = -2.0

    logger.debug(aligner)

    logger.info("> Aligning reference and consensus sequences")
    alignments = aligner.align(str(reference_record.seq), str(consensus_record.seq))
    alignment = alignments[0]

    target_locations = alignment.aligned[0] # Reference locations
    query_locations = alignment.aligned[1]  # Consensus locations

    logger.debug("ALIGNMENT: %s", alignment.aligned)

    with open("alignment.txt", "w") as f:
        f.write(str(alignment))

    # Account for the gaps in the alignment, by updating the consensus indexes
    # Needs to be double checked if the object is a tuple.
    logger.info("> Finding target tuples")
    indexes_differences = find_target_tuples_sorted(regions, target_locations)

    logger.info("> Updating query tuples")
    query_regions = update_query_tuple_with_difference(query_locations, indexes_differences)

    logger.info("> Updating consensus sequence")
    ref_seq = np.array(list(str(reference_record.seq)))
    cons_seq = np.array(list(str(consensus_record.seq)))

    cons_seq[query_regions] = ref_seq[regions]

    return "".join(cons_seq)


def find_target_tuples_sorted(regions, target_matrix):
    """
    Finds the target tuples sorted based on the regions.

    Args:
        regions (list): A sorted list of regions.
        target_matrix (matrix): A n x 2 matrix of [[start1,end1], [start2,end2], ...] alignment positions of the target query

    Returns:
        list: A list of tuples containing the index of the target tuple and the difference between the region and the start of the target tuple.
    """
    result = []

    for block_index, (start, end) in enumerate(target_matrix):
        logger.debug("target_matrix[%s]: start: %s - end: %s", block_index, start, end)

        # Iterate through regions while they are less than or equal to the end of the target tuple
        for region in regions:
            if region > end:
                logger.error("Region %s is greater than the end of the target tuple %s", region, end)
                logger.error("Please report this issue to the developers with your data from the current workdir.")
                break
            elif start <= region <= end:
                # If the region falls within the start and end of the target tuple,
                # Store the alignment block's index & the difference between the target position and the start alignment block.
                difference = region - start
                logger.debug("In alignment block (0-based) %s at postion %s (of that block), the consensus should be updated.", block_index, difference)
                result.append((block_index, difference))
    return result


def update_query_tuple_with_difference(query_matrix, results):
    updated_query_tuples = []

    for block_index, difference in results:
        if block_index is not None:
            query_tuple = query_matrix[block_index]
            updated_query_tuple = query_tuple[0] + difference
            updated_query_tuples.append(updated_query_tuple)

    return updated_query_tuples


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    reference = None
    consensus = None
    mpileup = None

    # Read in the reference sequence
    with open (args.reference, 'r') as f:
        reference = SeqIO.read(f, "fasta")
        logger.info("Reading reference ...\n")
    # Read in the consensus sequence
    with open (args.consensus, 'r') as f:
        consensus = SeqIO.read(f, "fasta")
        logger.info("Reading consensus ...\n")

    # Read in the mpileup file in a numpy array, Important to set the comments to None as '#' is used in the mpileup file
    with open (args.mpileup, 'r') as f:
        mpileup = np.loadtxt(f, dtype=str, delimiter="\t", comments=None)
        logger.info("Reading mpileup ...\n")

    # Check if mpileup is empty, if empty then exit
    if mpileup.size == 0:
        logger.error("Mpileup file is empty. Exiting ...")
        sys.exit(4)

    # Extract regions with coverage & subtract 1 for 0 index base
    low_coverage = mpileup[mpileup[:, 3].astype(int) <= args.minimum_depth, 1].astype(int) - 1

    # Annotate the consensus sequence with the reference sequence defined by low_coverage positions
    if len(low_coverage) > 0:
        logger.info("Low coverage regions found.")
        consensus_hybrid = annotate_ambiguous(reference, consensus, low_coverage, args)
    else:
        logger.info("No low coverage regions found. Writing out consensus sequence as is.")
        consensus_hybrid = consensus

    # Write out the annotated consensus sequence
    SeqIO.write(consensus_hybrid, f"{args.prefix}.fasta", "fasta")


if __name__ == "__main__":
    sys.exit(main())
