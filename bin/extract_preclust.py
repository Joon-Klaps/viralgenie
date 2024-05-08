#!/usr/bin/env python

import argparse
import logging
import sys

from Bio import SeqIO
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger()

taxon_map = {
        'Root':'R',
        'superkingdom':'D',
        'kingdom':'K',
        'phylum':'P',
        'class':'C',
        'order':'O',
        'family':'F',
        'genus':'G',
        'species':'S'}

TAXON_RANKED = {'R':9, 'D':8 , 'K':7, 'P':6, 'C':5, 'O':4, 'F':3, 'G':2, 'S':1}


def write_json(groups, prefix):
    """
    Write the groups to a json file.

    Args:
        groups (dict): A dictionary where the keys are taxids and the values are lists of rankedTaxon names.
        prefix (str): The prefix of the output file.
    """
    with open(f"{prefix}.json", "w", encoding="utf-8") as f_out:
        # Construct the JSON string manually
        json_str = "{\n"
        json_str += f'\t"ntaxa": {len(groups)},\n'
        json_str = json_str.rstrip(",\n") + "\n}"
        f_out.write(json_str)

def write_sequences(groups, sequences, prefix):
    """
    Extract the sequences from the input file based on the groups.

    Args:
        groups (dict): A dictionary where the keys are taxids and the values are lists of rankedTaxon objects.
        sequences (str): The path to the input sequence file.
        prefix (str): The prefix of the output file.
    """
    # check if groups has anything:
    sequence_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
    for taxid, ranked_taxon_list in groups.items():
        with open(f"{prefix}_taxid{taxid}.fa", "w", encoding='utf-8') as f_out:
            for ranked_taxon in ranked_taxon_list:
                if ranked_taxon.name in sequence_dict:
                    SeqIO.write(sequence_dict[ranked_taxon.name], f_out, "fasta")
                else:
                    logger.warning("The sequence %s was not found in the input sequence file %s!", ranked_taxon.name, sequences)

def write_report(groups, prefix):
    """
    Write the groups to a report file.
    """

    with open(f"{prefix}.resolved.txt", "w", encoding="utf-8") as f_out:

        for taxid, ranked_taxon_list in groups.items():
            for ranked_taxon in ranked_taxon_list:
                f_out.write(f"{ranked_taxon}\n")

def simplify_taxonomic_ranks(dic, nodes, rank):
    """
    Simplify the taxonomic ranks of the given dictionary up to a certain level.

    Args:
        dic: The dictionary of taxids and a list of RankedTaxon objects.
        nodes: The dictionary of taxonomic nodes.
        rank: The rank to simplify the taxonomic ranks up to.

    Returns:
        Dict: A dictionary with taxid and a list containing RankedTaxon objects referring to the resolved merged reads of that taxid.
    """

    results_dic = defaultdict(list)
    for taxid, ranked_taxon_list in dic.items():
        node = nodes[taxid]
        result = node.get_parent_of_rank(rank)

        if taxid != result.taxid:
            logger.debug("Simplifying taxonomic ranks of taxid %s to %s", taxid, result.taxid)
            # Create a copy of the ranked_taxon_list to avoid modifying the original
            simplified_list = list(ranked_taxon_list)
            # update all the taxon ids in the copy
            for ranked_taxon in simplified_list:
                ranked_taxon.taxid = result.taxid

            results_dic[result.taxid].extend(simplified_list)

        else:
            # No change in taxa, simply copy the list (optional)
            results_dic[result.taxid].extend(list(ranked_taxon_list))

    logger.info("Classifications resolved successfully of %d hit(s)!", get_group_size(results_dic))
    return results_dic

def file_not_found(file):
    """
    Log an error message and exit the program if a file is not found.
    """
    logger.error("The given input file %s, was not found!", file)
    sys.exit(2)

def first_common_element(x,y):
    """
    Find the first common element between two lists.

    Args:
        x: The first list.
        y: The second list.

    Returns:
        The first common element between the two lists.
    """
    for i in x:
        if i in y:
            return i
    return None

def resolve_read_classifications(args, nodes):
    """
    resolve the read classifications from Kaiju and Kraken.

    Args:
        args: The parsed command line arguments.
        nodes: The dictionary of taxonomic nodes.

    Returns:
        Dict: A dictionary with taxid and a list containing RankedTaxon objects referring to the resolved merged reads of that taxid.
    """

    resolved_results = defaultdict(list)

    logger.info("Resolving the classifications ...")

    # Use Kaiju and Kraken
    if args.kaiju_classifications and args.kraken_classifications:

        # Open files and iterate over every line
        with open (args.kaiju_classifications, encoding="utf-8") as kaiju_file, open(args.kraken_classifications, encoding="utf-8") as kraken_file:
            for line_number, (kaiju_line, kraken_line) in enumerate(zip(kaiju_file, kraken_file)):
                kaiju_values =kaiju_line.strip().split("\t")[0:3]
                kraken_values = kraken_line.strip().split("\t")[0:3]

                # logger.debug("Looking at Kaiju - Kraken read classifications (%s - %s) at line %d", kaiju_values, kraken_values, line_number)
                if len(kaiju_values) != 3 or len(kraken_values) != 3:
                    logger.error("The Kaiju - Kraken read classifications (%s - %s) do not have the same number of columns at line %d", kaiju_values, kraken_values, line_number)
                    sys.exit(1)

                # Define a new conflict object
                conflict = ConflictTaxa(kaiju_values = kaiju_values, kraken_values = kraken_values)

                # Check if it's about the same read
                if conflict.kaiju.name != conflict.kraken.name:
                    logger.error (
                        "The Kaiju - Kraken read names (%s - %s)  do not match at line %d",
                        conflict.kaiju.name,
                        conflict.kraken.name,
                        line_number
                    )
                    logger.error (
                        "PLEASE RESUBMIT with --kaiju-classifications <(sort -k2,2 <kaiju-file>) --kraken-classifications <(sort -k2,2 <kraken-file>)"
                    )
                    sys.exit(1)

                # Resolve the read using the specified mere strategy
                if args.merge_strategy == "1":
                    conflict.use_kaiju()
                elif args.merge_strategy == "2":
                    conflict.use_kraken()
                elif args.merge_strategy == "lca":
                    conflict.use_lca(nodes)
                elif args.merge_strategy == "lowest":
                    conflict.use_lowest(nodes)
                else :
                    logger.error("Unknown merge strategy specified %s", args.merge_strategy )

                # Store result
                if not conflict.result:
                    logger.error("Something went wrong, no result found for the conflict at line %d", line_number)
                    sys.exit(1)
                else: # append the result to the resolved results
                    # logger.debug("Resolved read classification %s for line %d", conflict.result, line_number)
                    resolved_results[conflict.result.taxid].append(conflict.result)

    # Use Kaiju
    elif args.kaiju_classifications:
        with open (args.kaiju_classifications, encoding="utf-8") as kaiju_file:
            for line_number, kaiju_line in enumerate(kaiju_file):
                kaiju_values =kaiju_line.strip().split("\t")[0:3]

                logger.debug("Looking at Kaiju read classifications (%s) at line %d", kaiju_values, line_number)
                if len(kaiju_values) != 3:
                    logger.error("The Kaiju read classifications (%s) do not have the same number of columns at line %d", kaiju_values, line_number)
                    sys.exit(1)

                # Define a new conflict object
                conflict = ConflictTaxa(kaiju_values = kaiju_values, kraken_values = None)
                conflict.use_kaiju()

                if not conflict.result:
                    logger.error("Something went wrong, no result found for the conflict at line %d", line_number)
                    sys.exit(1)
                else: # append the result to the resolved results
                    logger.debug("Resolved read classification %s for line %d", conflict.result, line_number)
                    resolved_results[conflict.result.taxid].append(conflict.result)

    # Use Kraken
    else:
        with open(args.kraken_classifications, encoding="utf-8") as kraken_file:
            for line_number, kraken_line in enumerate(kraken_file):
                kraken_values = kraken_line.strip().split("\t")[0:3]

                logger.debug("Looking at Kraken read classifications (%s) at line %d", kraken_values, line_number)
                if len(kraken_values) != 3:
                    logger.error("The Kraken read classifications (%s) do not have the same number of columns at line %d", kraken_values, line_number)
                    sys.exit(1)

                # Define a new conflict object
                conflict = ConflictTaxa(kaiju_values = None, kraken_values = kraken_values)
                conflict.use_kraken()

                if not conflict.result:
                    logger.error("Something went wrong, no result found for the conflict at line %d", line_number)
                    sys.exit(1)
                else: # append the result to the resolved results
                    logger.debug("Resolved read classification %s for line %d", conflict.result, line_number)
                    resolved_results[conflict.result.taxid].append(conflict.result)

    logger.info("Classifications resolved successfully of %d hit(s)!", get_group_size(resolved_results))
    return resolved_results

def deduplicate(x):
    """
    Remove duplicates from a list while preserving the order.
    """
    if x:
        logger.debug(x)
        return list(dict.fromkeys(tuple(x)))
    return x

def define_lists(args, nodes):
    """
    Define the whitelist and blacklist based on the command line arguments include|exclude children|parents.

    Args:
        args: The parsed command line arguments.
        nodes: The dictionary of taxonomic nodes.

    Returns:
        tuple: A tuple containing the whitelist and blacklist.
    """
    whitelist = []
    blacklist = []

    # whitelist
    if args.include_children:
        whitelist.extend(args.include_children)
        for taxid in args.include_children:
            if taxid not in nodes:
                logger.warning("The taxid %s was not found in the taxonomy file!", taxid)
            else:
                whitelist.extend(nodes[taxid].get_all_child_taxids())
    if args.include_parents:
        whitelist.extend(args.include_parents)
        for taxid in args.include_parents:
            if taxid not in nodes:
                logger.warning("The taxid %s was not found in the taxonomy file!", taxid)
            else:
                whitelist.extend(nodes[taxid].get_all_parent_taxids())

    # blacklist
    if args.exclude_children:
        blacklist.extend(args.exclude_children)
        for taxid in args.exclude_children:
            if taxid not in nodes:
                logger.warning("The taxid %s was not found in the taxonomy file!", taxid)
            else:
                blacklist.extend(nodes[taxid].get_all_child_taxids())
    if args.exclude_parents:
        blacklist.extend(args.exclude_parents)
        for taxid in args.exclude_parents:
            if taxid not in nodes:
                logger.warning("The taxid %s was not found in the taxonomy file!", taxid)
            else:
                blacklist.extend(nodes[taxid].get_all_parent_taxids())

    whitelist = deduplicate(whitelist)
    blacklist = deduplicate(blacklist)

    overlap = set(whitelist) & set(blacklist)
    if overlap:
        logger.warning("Overlap between 'to include taxa' and 'to remove taxa', removing them from 'to remove list': %s \n\t TIP: If you want to really remove this taxa, specify '--exclude-taxa' ", overlap)
        blacklist = list(set(blacklist).difference(overlap))

    if args.exclude_taxa:
        blacklist.extend(args.exclude_taxa)
        blacklist = deduplicate(blacklist)

    if whitelist:
        logger.info("White list of taxids created, size: %d", len(whitelist))

    if blacklist:
        logger.info("Black list of taxids created, size: %d", len(blacklist))

    return whitelist, blacklist

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
        if level_rank not in taxon_map:
            level_rank = '-'
        else:
            level_rank = taxon_map[level_rank]
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

    logger.warning("Using kraken report to reconstruct taxonomy tree, consider providing a nodes.dmp file with --nodes for complete taxonomy tree.")

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
                root_node = Tree(taxid = str(taxid), level_rank= level_rank[0], level_num= level_num)
                prev_node = root_node
                continue

            #move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent

            #determine correct level ID
            if level_rank == '-' or len(level_rank) > 1:
                if prev_node.level_rank in taxon_map.values():
                    level_rank = prev_node.level_rank + '1'
                else:
                    num = int(prev_node.level_rank[-1]) + 1
                    level_rank = prev_node.level_rank[:-1] + str(num)

            curr_node = Tree(taxid = str(taxid), level_rank = level_rank[0], level_num=level_num, parent=prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node
            nodes[str(taxid)] = curr_node

    return nodes

def get_group_size(group):
    """
    Get the size of a group of reads.

    Args:
        group: A list of reads.

    Returns:
        The size of the group.
    """
    size = 0
    for key, taxon_list in group.items():
        size += len(taxon_list)
    return size

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

    nodes['0'] = Tree(taxid = '0', level_rank = 'R')

    with open(nodes_file, encoding= 'utf-8') as f:
        logger.info("0 nodes read")
        for line in f:
            count_nodes += 1
            if (count_nodes % 500000) == 0:
                logger.info("%d nodes read", count_nodes)

            [curr_taxid,parent_taxid,rank] = line.strip().split('\t|\t')[0:3]

            new_rank = '-'
            if rank in taxon_map:
                new_rank = taxon_map[rank]

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

def process_taxonomy(nodes_file, kraken_report):
    """
    Process the taxonomy file based on the command line arguments.

    Args:
        nodes_file: the parsed NCBI nodes file location
        kraken_report: the parsed Kraken report file location

    Returns:
        dict: A dictionary where the keys are tax_ids and values are dictionaries containing information like parent_tax_id and rank (if provided).
    """
    if nodes_file:
        if not nodes_file.is_file():
            file_not_found(nodes_file)
        nodes = parse_nodes_dmp(nodes_file)
    elif kraken_report:
        if not kraken_report.is_file():
            file_not_found(kraken_report)
        nodes = parse_kraken_report(kraken_report)
    else:
        logger.error("Please provide either an '--nodes' or '--database' or '--kraken_report' as %s or %s was not found!", nodes, kraken_report)
        sys.exit(2)
    logger.info("Taxonomy read in with %d taxa", len(nodes))
    return nodes

# Taken and modified from https://github.com/jenniferlu717/KrakenTools/blob/master/make_ktaxonomy.py
class Tree:
    """
    A class to represent a tree of taxonomic ranks. Containing information on parent and their children
    """

    def __init__(self, taxid, level_rank, level_num=-1, parent=None, children=None):
        self.taxid = taxid
        self.level_rank = level_rank
        # Other attributes for later
        self.name = ""
        self.level_num = level_num
        self.p_taxid = -1
        # Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        """
        Extend the current child node list to the current node
        """
        assert isinstance(node, Tree)
        self.children.append(node)

    def _get_taxid(self):
        """
        Return the taxid of the current node
        """
        return self.taxid

    def __str__(self):
        parent_taxid = None
        if self.parent:  # Check if parent exists before accessing its attribute
            parent_taxid = self.parent._get_taxid()
        return f"taxid:{self.taxid} - rank:{self.level_rank} - parent:{parent_taxid} - number of children: {len(self.children)}"

    def get_all_child_taxids(self):
        """
        Traverses the tree structure recursively, collecting taxids of all child nodes.

        Args:
            self: The current node of the tree structure.

        Returns:
            A list containing all child taxids encountered during traversal.
        """

        all_taxids = []

        def traverse_children(node):
            if not node.children:
                return  # Base case: Leaf node (no children)

            for child in node.children:
                all_taxids.append(child._get_taxid())  # Add child's taxid
                traverse_children(child)  # Recursive call for each child

        traverse_children(self)
        return all_taxids

    def get_all_parent_taxids(self):
        """
        Traverses the tree structure recursively, collecting taxids of all parent nodes.

        Args:
            self: The current node of the tree structure.

        Returns:
            A list containing all parent taxids encountered during traversal.
        """

        all_taxids = []

        def traverse_parents(node):
            all_taxids.append(node._get_taxid())
            if not node.parent:
                return
            traverse_parents(node.parent)

        traverse_parents(self)
        return all_taxids

    def get_parent_of_rank(self, rank):
        """
        Traverses the tree structure recursively, collecting taxids of all parent nodes.

        Args:
            self: The current node of the tree structure.
            rank: The rank of the parent node to search for.

        Returns:
            The parent Tree object that holds the specified rank or self if the rank is not found or if self is higher then rank.
        """
        self_level = TAXON_RANKED.get(self.level_rank, -1)
        rank_level = TAXON_RANKED.get(rank, -1) #Should never return -1 though

        if rank_level == -1:
            logger.error("The rank %s is not a valid rank", rank)
            sys.exit(1)

        if rank_level <= self_level:
            return self

        # Should normally never happen that we simplify up to root
        if self.parent is None:
            logger.warning("Reached non-existend parent of %s", self.taxid)
            return self  # Reached the root node, no higher parent

        return self.parent.get_parent_of_rank(rank)

class RankedTaxon:
    """
    A class to represent the first three elements of a Kaiju or Kraken hit.
    """
    def __init__(self, values):
        logger.debug("Parsing values %s", values)
        self.classified = values[0]
        self.name = values[1]
        self.taxid = values[2]

    def __str__(self):
        return f"{self.classified}\t{self.name}\t{self.taxid}"
class ConflictTaxa:
    """
    A class to represent the combination of Kaiju & Kraken and their result.
    """
    def __init__(self, kaiju_values, kraken_values):
        self.kaiju  = None
        self.kraken = None
        self.result = None

        if kaiju_values:
            self.kaiju = RankedTaxon(kaiju_values)
        if kraken_values:
            self.kraken = RankedTaxon(kraken_values)

    def _has_conflict(self):
        """
        Check if the Kaiju and Kraken classifications are different.
        """
        if self.kaiju.classified == "U" and self.kraken.classified == "C":
            self.result = self.kraken
            return False
        elif self.kaiju.classified == "C" and self.kraken.classified == "U":
            self.result = self.kaiju
            return False
        elif self.kaiju.taxid == self.kraken.taxid:
            self.result = self.kaiju
            return False
        return True


    def use_kaiju(self):
        """
        Sets the result to the Kaiju classification.
        """
        # make it so we can use the kaiju classification as the result without checking the kraken classification
        if not isinstance(self.kraken, RankedTaxon):
            self.result = self.kaiju
        elif self.kaiju.classified == "U" and self.kraken.classified == "C":
            self.result = self.kraken
        else:
            self.result = self.kaiju

    def use_kraken(self):
        """
        Sets the result to the Kraken classification.
        """
        # make it so we can use the kraken classification as the result without checking the kaiju classification
        if not isinstance(self.kaiju, RankedTaxon):
            self.result = self.kraken
        elif self.kaiju.classified == "C" and self.kraken.classified == "U":
            self.result = self.kaiju
        else:
            self.result = self.kraken


    def use_lowest(self, nodes):
        """
        Sets the result to the lowest ranking of the two taxon identifiers.
        """
        if not self._has_conflict():
            return

        # check if kaiju parent of Kraken
        if self.kaiju.taxid in nodes[self.kraken.taxid].get_all_parent_taxids():
            logger.debug("Kraken is lower: %s - %s", nodes[self.kraken.taxid] , nodes[self.kaiju.taxid])
            self.result = self.kraken
            return

        # check if kraken parent of Kaiju
        if self.kraken.taxid in nodes[self.kaiju.taxid].get_all_parent_taxids():
            logger.debug("Kaiju is lower: %s - %s", nodes[self.kaiju.taxid] , nodes[self.kraken.taxid])
            self.result = self.kaiju
            return

        # if one is not a parent of the other, use LCA
        self.use_lca(nodes)


    def use_lca(self, nodes):
        """
        Sets the result to the lowest common ancestor of the two taxon identifiers.
        """
        if not self._has_conflict():
            return
        kaiju_parents = nodes[self.kaiju.taxid].get_all_parent_taxids()
        kraken_parents = nodes[self.kraken.taxid].get_all_parent_taxids()

        common_parent = first_common_element(kaiju_parents, kraken_parents)
        if common_parent:
            result = RankedTaxon(["C", self.kaiju.name, common_parent])
            self.result = result
        else:
            logger.error("No common parent found between Kaiju and Kraken classifications: %s - %s", self.kaiju, self.kraken)
            sys.exit(1)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract the sequences based on the results of Kaiju, Kraken or both of them combined.",
        epilog="Example: python extract_precluster.py input.tsv sequence.fa prefix",
    )

    parser.add_argument(
        "--sequences",
        metavar="SEQUENCES",
        required=True,
        type=Path,
        help="Input sequence file used for looking up",
    )

    parser.add_argument(
        "-p",
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
        "-db",
        "--database",
        metavar="DATABASE",
        dest="database",
        type=Path,
        help="A database directory containing an NCBI nodes file. `nodes.dmp`",
    )

    parser.add_argument(
        "--kraken-report",
        metavar="KRAKEN_REPORT",
        type=Path,
        help="Kraken's report `--report` containing information on the the parent taxa",
    )

    parser.add_argument(
        "-c",
        "--merge-strategy",
        metavar="MERGE-STRATEGY",
        type=str,
        dest="merge_strategy",
        help="Specify the merge strategy of the classifications of Kaiju and Kraken. '1' for Kaiju, '2' for Kraken, 'lca' for lowest common ancestor and 'lowest' for lowest ranking of the two taxon identifiers.",
        choices=("1", "2", "lca", "lowest"),
        default="lca",
    )

    parser.add_argument(
        "-ic",
        "--include-children",
        nargs="+",
        type=str,
        dest = "include_children",
        help="A list of taxids to whitelist during filtering and include their children.",
    )

    parser.add_argument(
        "-ec",
        "--exclude-children",
        nargs="+",
        type=str,
        dest = "exclude_children",
        help="A list of taxids to blacklist during filtering and exclude their children.",
    )

    parser.add_argument(
        "-ip",
        "--include-parents",
        nargs="+",
        type=str,
        dest = "include_parents",
        help="A list of taxids to whitelist during filtering and include their parents.",
    )

    parser.add_argument(
        "-ep",
        "--exclude-parents",
        nargs="+",
        type=str,
        dest = "exclude_parents",
        help="A list of taxids to blacklist during filtering and exclude their parents.",
    )

    parser.add_argument(
        "-r",
        "--exclude-taxa",
        nargs="+",
        type=str,
        dest = "exclude_taxa",
        help="A list of taxids to blacklist during filtering.",
    )

    parser.add_argument(
        "-s",
        "--simplification-level",
        type= str,
        dest = "simplification_level",
        help="The level of simplification of the taxonomic ranks. 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'.",
        choices=("superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"),
    )

    parser.add_argument(
        "-u",
        "--keep-unclassified",
        action="store_true",
        default= False,
        help="Keep unclassified reads in the output.",
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

    if not args.kaiju_classifications and not args.kraken_classifications:
        logger.error("Please provide either an '--kaiju_classifications' or '--kraken_classifications' as %s or %s was not found!", args.kaiju_classifications, args.kraken_classifications)
        sys.exit(2)

    need_taxonomy = (
        args.simplification_level or
        args.exclude_children or
        args.include_children or
        args.exclude_parents or
        args.include_parents or
        args.merge_strategy in ['lowest', 'lca']
    )
    nodes = {}
    if need_taxonomy:
        logger.info("Processing taxonomy file")
        if args.database:
            args.nodes = args.database / "nodes.dmp"
        nodes = process_taxonomy(args.nodes, args.kraken_report)

    # A dictionary of taxid and a list containing RankedTaxon objects referring to the resolved merged reads of that taxid.
    resolved = resolve_read_classifications(args, nodes)
    results = resolved.copy()

    # Keeping only the specified taxids
    whitelist, blacklist = define_lists(args, nodes)
    if whitelist:
        logger.info("Keeping only whitlisted taxids")
        results = {key: value for key, value in results.items() if key in whitelist}
        logger.info("Kept only %d taxids", get_group_size(results))

    if blacklist:
        logger.info("Removing blacklisted taxids")
        tmp = get_group_size(results)
        results = {key: value for key, value in results.items() if key not in blacklist}
        logger.info("Removed %d hits", (tmp-get_group_size(results)))

    if not args.keep_unclassified:
        logger.info("Removing unclassified hits")
        tmp = get_group_size(results)
        results = {key: value for key, value in results.items() if key != "0"}
        logger.info("Removed %d hits", (tmp-get_group_size(results)))

    # Reduce the taxonomic ranks up to a certain level
    if args.simplification_level:
        logger.info("Simplifying taxonomic ranks to %s", args.simplification_level)
        tmp  = len(results)
        results = simplify_taxonomic_ranks(results, nodes, taxon_map[args.simplification_level])
        logger.info("Simplified from %d to %d taxonomie(s)", tmp, len(results))

    if not results:
        logger.error("No final taxon groups found, is it possible everything was removed? exiting!")
        logger.error("Resolved taxons: %s ", resolved.keys())
        if whitelist:
            logger.error("Whitelist (keep only) taxons: %s", whitelist)
        if blacklist:
            logger.error("Blacklist (remove) taxons: %s", blacklist)
        sys.exit(1)

    logger.info("Writing sequences and json file")
    write_sequences(results, args.sequences, args.prefix)
    write_json(results, args.prefix)
    write_report(results, args.prefix)

    logger.info("All done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
