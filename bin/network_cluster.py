#!/usr/bin/env python

"""Cluster based on a the created network generated from tools"""

import argparse
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import re
import logging
import sys
from pathlib import Path

logger = logging.getLogger()

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to extract sequence names from cdhit's cluster files.",
        epilog="Example: python blast_filter.py in.clstr prefix",
    )

    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="cluster file from chdit or vsearch containing cluster information.",
    )

    parser.add_argument(
        "method",
        metavar="METHOD",
        type=Path,
        help="Comparison method containing the necessary information for creating a network.",
    )

    parser.add_argument(
        "file_out_prefix",
        metavar="FILE_OUT_PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "-s",
        "--score",
        metavar="score",
        type=float,
        help="Score cutoff for clustering",
        default=0.80,
    )

    parser.add_argument(
        "-p",
        "--pattern",
        metavar="PATTERN",
        type=str,
        help="Regex pattern to filter clusters by centroid sequence name.",
        default="^(TRINITY)|(NODE)|(k\d+)",  # Default pattern matches Trinity, SPADes and MEGAHIT assembly names
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)

def read_in_file(file_in, method):
    """
    Read in the file and return a networkx graph object
    """
    method_dict = {
        'mash': read_in_mash
        # Add more methods here if needed
    }

    # Check if the method exists in the dictionary
    if method in method_dict:
        return method_dict[method](file_in)
    else:
        raise ValueError(f"Method '{method}' not found.")

def read_in_mash(file_in):
    """
    Read in the file and return a networkx graph object
    """
    df = pd.read_csv(file_in,sep="\t", index_col="#query")
    G = nx.from_pandas_adjacency(df)
    return G


def filter_network(network, treshold):
    """
    Filter the network based on the given score
    """
    filtered_network = network.copy()
    edges_to_remove = [(u, v) for u, v, d in filtered_network.edges(data=True) if d["weight"] >= threshold ]
    return filtered_network.remove_edges_from(edges_to_remove)

def cluster_network(network):
    """
    Cluster the network based on the given score
    """
    partition = nx.community.louvain_communities(network, seed=10)
    return partition

def to_tsv(network, file_out_prefix):
    """
    Write the network to a tsv file
    """
    # Create a list of tuples containing the first word (if it's a set) and its index
    indexed_data = [(list(s)[0], i) for i, s in enumerate(data) if s]

    # Write the indexed data to a TSV file
    with open(f'{file_out_prefix}.tsv', 'w') as file:
        for line in indexed_data:
            file.write(f"{line[0]}\t{line[1]}\n")

# have a look at this example on how to do this.
# https://community.plotly.com/t/displaying-edge-labels-of-networkx-graph-in-plotly/39113/2
# Also go into more depth on how to visualize with weights in networkx or use an heatmap
#def visualize_network(network):

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)

    network = read_in_file(args.file_in, args.method)

    network = filter_network(network, 1 - args.score)

    clusters = cluster_network(network)

    to_tsv(clusters, args.file_out_prefix)

if __name__ == "__main__":
    sys.exit(main())
