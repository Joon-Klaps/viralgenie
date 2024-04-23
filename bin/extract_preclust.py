#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path

# from Bio import SeqIO

logger = logging.getLogger()

map_ranks = {
        'Root':'R',
        'kingdom':'K',
        'superkingdom':'D',
        'phylum':'P',
        'class':'C',
        'order':'O',
        'family':'F',
        'genus':'G',
        'species':'S'}

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

def process_kraken_report(report_line):
    """
    Parses single line from report output and returns taxID, levelID
    Taken from https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py

    Args:
        kraken report file with the following tab delimited lines
            - percent of total reads
            - number of reads (including at lower levels)
            - number of reads (only at this level)
            - taxonomy classification of level
                (U, - (root), - (cellular org), D, P, C, O, F, G, S)
            - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
            - spaces + name
    Returns:
        A list containing the following elements:
            - taxonomy ID
            - level number (number of spaces before name)
            - level_rank (type of taxonomy level - U, R, D, P, C, O, F, G, S, etc)
    """
    l_vals = report_line.strip().split('\t')
    if len(l_vals) < 5:
        return []
    try:
        int(l_vals[1])
    except ValueError:
        return []
    #Extract relevant information
    try:
        taxid = int(l_vals[-3])
        level_rank = l_vals[-2]
        if level_rank not in map_ranks:
            level_rank = '-'
        else:
            level_rank = map_ranks[level_rank]
    except ValueError:
        taxid = int(l_vals[-2])
        level_rank = l_vals[-3]
    #Get spaces to determine level num
    spaces = 0
    for char in l_vals[-1]:
        if char == ' ':
            spaces += 1
        else:
            break
    level_num = int(spaces/2)
    return[taxid, level_num, level_rank]

def parse_kraken_report(kraken_report):
    """
    Parses a Kraken report file and returns a dictionary of nodes.
    Based on https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py

    Args:
        kraken_report: A file object containing the Kraken report.

    Returns:
        A dictionary where keys are tax_ids and values are dictionaries containing information like parent_tax_id and rank (if provided).
    """
    nodes = {}
    with open(kraken_report, encoding="utf-8") as f:
        prev_node = -1
        for line in f:
            #extract values
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            if len(report_vals) != 3:
                logger.error("Corrupted kraken report line: %s", line)
            taxid, level_num, *level_rank = report_vals

            if taxid == 0:
                continue

            if taxid == 1:
                level_rank = 'R'
                root_node = Tree(taxid = taxid, level_rank= level_rank, level_num= level_num)
                prev_node = root_node
                continue

            #move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent

            #determine correct level ID
            if level_rank == '-' or len(level_rank) > 1:
                if prev_node.level_rank in map_ranks.values():
                    level_rank = prev_node.level_rank + '1'
                else:
                    num = int(prev_node.level_rank[-1]) + 1
                    level_rank = prev_node.level_rank[:-1] + str(num)

            curr_node = Tree(taxid = taxid, level_rank = level_rank, level_num=level_num, parent=prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node
            nodes[taxid] = curr_node

    return nodes

def parse_nodes_dmp(nodes_file):
    """
    Parses a file containing variable-format node information and returns two dictionaries.

    Args:
        nodes_file: An open file object in text mode containing node data.

    Returns:
        A dictionary where keys are tax_ids and values are dictionaries containing information like parent_tax_id and rank (if provided).
    """
    nodes = {}
    p_notsaved = {}
    count_nodes = 0
    with open(nodes_file, encoding= 'utf-8') as f:
        logger.debug("0 nodes read")
        for line in f:
            count_nodes += 1
            if (count_nodes % 500000) == 0:
                logger.debug("%d nodes read", count_nodes)

            [curr_taxid,parent_taxid,rank] = line.strip().split('\t|\t')[0:3]

            new_rank = '-'
            if rank in map_ranks:
                new_rank = map_ranks[rank]

            curr_node = Tree(taxid = curr_taxid, level_rank = new_rank)
            nodes[curr_taxid] = curr_node

            if curr_taxid == "1":
                curr_node.level_rank = 'R'

            elif parent_taxid in nodes:
                #save parent
                curr_node.parent = nodes[parent_taxid]
                curr_node.p_taxid = parent_taxid
                nodes[parent_taxid].add_child(curr_node)
            else:
                #parent not linked
                p_notsaved[curr_taxid] = curr_node
                curr_node.p_taxid = parent_taxid

    for taxid, node  in p_notsaved.items():
        p_taxid = node.p_taxid
        if p_taxid not in nodes:
            logger.error("ERROR corrupted nodes.dmp: %s not found in nodes.dmp file", p_taxid)
            continue
        p_node = nodes[p_taxid]
        node.parent = p_node
        p_node.add_child(node)

    return nodes


# Taken and modified from https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
class Tree:
    """
    A class to represent a tree of taxonomic ranks. Containing information on parent and their children
    """
    def __init__(self,  taxid, level_rank, level_num=-1, parent=None,children=None):
        self.taxid = taxid
        self.level_rank= level_rank
        #Other attributes for later
        self.name = ''
        self.level_num = level_num
        self.p_taxid = -1
        #Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        """
        Extend the current child node list to the current node
        """
        assert isinstance(node,Tree)
        self.children.append(node)
    def _get_taxid(self):
        """
        Return the taxid of the current node
        """
        return self.taxid

    def __str__(self):
        return f"taxid:{self.taxid} - rank:{self.level_rank} - parent:{self.parent._get_taxid()}"

class RankedTaxon:
    """
    A class to represent a ranked taxon of Kaiju and (/or) Kraken hit. Containing information on parent
    """
    def __init__(self, kaiju_taxid, kraken_taxid, classified, name):
        self.kaiju_taxid = kaiju_taxid
        self.kraken_taxid = kraken_taxid
        self.classified = classified
        self.name = name
        self.taxon_rank = -1
        self.result_taxid = -1 # can either be the result from the merging or the simplified version after where it is collapsed up to a specific taxon level.

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
        "--sequences",
        metavar="SEQUENCES",
        # required=True,
        type=Path,
        help="Input sequence file used for looking up",
    )

    parser.add_argument(
        "-p"
        "--prefix",
        metavar="PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "--kaiju-classifications",
        metavar="KAIJU_CLASSIFICATIONS",
        type=Path,
        help="Classified contigs|reads using Kaiju with only the first 3 columns sorted by readID.",
    )

    parser.add_argument(
        "--kraken-classifications",
        metavar="KRAKEN_CLASSIFICATIONS",
        type=Path,
        help="Classified contigs|reads using Kraken with only the first 3 columns sorted by readID.",
    )

    parser.add_argument(
        "-n",
        "--nodes",
        metavar="NODES",
        type=Path,
        help=" NCBI nodes file. `nodes.dmp`",
    )

    parser.add_argument(
        "--kraken-report",
        metavar="KRAKEN_REPORT",
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
        default="INFO",
    )
    return parser.parse_args(argv)

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level,
        format="[%(asctime)s - %(levelname)s] %(message)s",
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    # if not args.kaiju_classifications.is_file() and not args.kraken_classifications.is_file():
    #     logger.error("Please provide either an '--kaiju_classifications' or '--kraken_classifications' as %s or %s was not found!", args.kaiju_classifications, args.kraken_classifications)
    #     sys.exit(2)
    # if not args.classifications.is_file():
    #     logger.error("The given input file %s was not found!", args.classifications)
    #     sys.exit(2)
    # if not args.sequences.is_file():
    #     logger.error("The given input file %s was not found!", args.sequences)
    #     sys.exit(2)

    if args.nodes:
        nodes = parse_nodes_dmp(args.nodes)
    if args.kraken_report:
        nodes = parse_kraken_report(args.kraken_report)


    counter = 0
    for key, value in nodes.items():
        counter+=1
        print(f"Key: {key}, Value: {value}")
        if counter == 5:  # Check if we've iterated through 5 elements based on the size of OrderedDict
            break

    # groups = create_groups(args.classifications)

    # extract_sequences(groups, args.sequences, args.file_out_prefix)
    # write_json(groups, args.file_out_prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
