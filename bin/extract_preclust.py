#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger()




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
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split("\t")
            classified = parts[0]
            seq_name = parts[1]
            taxid = parts[2]
            if classified == "C":
                if taxid not in groups:
                    groups[taxid] = []
                groups[taxid].append(seq_name)
            elif classified == "U":
                if "U" not in groups:
                    groups["U"] = []
                groups["U"].append(seq_name)
    return groups


def write_json(groups, file_out_prefix):
    """
    Write the groups to a json file.

    Args:
        groups (dict): A dictionary where the keys are taxids and the values are lists of sequence names.
            The 'U' key represents unclassified sequences.
        file_out_prefix (str): The prefix of the output file.
    """
    with open(f"{file_out_prefix}.json", "w") as f_out:
        # Construct the JSON string manually
        json_str = "{\n"
        json_str += f'\t"ntaxa": {len(groups)},\n'
        json_str = json_str.rstrip(",\n") + "\n}"
        f_out.write(json_str)


def extract_sequences(groups, sequences, file_out_prefix):
    """
    Extract the sequences from the input file based on the groups.

    Args:
        groups (dict): A dictionary where the keys are taxids and the values are lists of sequence names.
            The 'U' key represents unclassified sequences.
        sequences (str): The path to the input sequence file.
        file_out_prefix (str): The prefix of the output file.
    """
    sequence_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
    for taxid, seq_names in groups.items():
        with open(f"{file_out_prefix}_taxid{taxid}.fa", "w") as f_out:
            for seq_name in seq_names:
                if seq_name in sequence_dict:
                    SeqIO.write(sequence_dict[seq_name], f_out, "fasta")

class rankedTaxon:
    """
    A class to represent a ranked taxon of Kaiju and (/or) Kraken hit. Containing information on parent
    """
    def __init__(self, kaiju_taxid, kraken_taxid, classified, name):
        self.kaiju_taxid = kaiju_taxid
        self.kraken_taxid = kraken_taxid
        self.classified = classified
        self.name = name
        self.taxon_rank = self._lookup_taxon_rank()

    def _lookup_parent_taxa(self):
        return 0

    def _lookup_taxon_rank(self):
        return 0

    def _merge_classification(self, merge_strategy):
        return 0

    def _simplify_taxon_rank(self, simplification_level):
        return 0

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
        "kaiju_classifications",
        metavar="KAIJU_CLASSIFICATIONS",
        type=Path,
        help="Classified contigs|reads using Kaiju with only the first 3 columns sorted by readID.",
    )

    parser.add_argument(
        "kraken_classifications",
        metavar="KRAKEN_CLASSIFICATIONS",
        type=Path,
        help="Classified contigs|reads using Kraken with only the first 3 columns sorted by readID.",
    )

    parser.add_argument(
        "kaiju_nodes",
        metavar="KAIJU_NODES",
        type=Path,
        help="Kaiju nodes file. `nodes.dmp`",
    )

    parser.add_argument(
        "kraken_report",
        metavar="KAIJU_NAMES",
        type=Path,
        help="Kraken's report `--report` containing information on the the parent taxa",
    )

    parser.add_argument(
        "-c"
        "--merge-classifications",
        type=str,
        help="Specify the merge strategy of the classifications of Kaiju and Kraken. '1' for Kaiju, '2' for Kraken, 'lca' for lowest common ancestor and 'lowest' for lowest ranking of the two taxon identifiers.",
        choices=("1", "2", "lca", "lowest"),
        default="lca",
    )

    parser.add_argument(
        "-ic",
        "--include-children",
        nargs="+",
        type=str,
        help="A list of taxids to whitelist during filtering and include their children.",
    )

    parser.add_argument(
        "-ec",
        "--exclude-children",
        nargs="+",
        type=str,
        help="A list of taxids to blacklist during filtering and exclude their children.",
    )

    parser.add_argument(
        "-ip",
        "--include-parents",
        nargs="+",
        type=str,
        help="A list of taxids to whitelist during filtering and include their parents.",
    )

    parser.add_argument(
        "-ep",
        "--exclude-parents",
        nargs="+",
        type=str,
        help="A list of taxids to blacklist during filtering and exclude their parents.",
    )

    parser.add_argument(
        "-s",
        "--simplification-level",
        type= str,
        help="The level of simplification of the taxonomic ranks. 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'.",
        choices=("superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"),
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)

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
    write_json(groups, args.file_out_prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
