#!/usr/bin/env python

"""Cluster based on a the created network generated from tools"""

import argparse
import igraph as ig
import pandas as pd
import leidenalg as la
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to create clusters or communities from distance measures using the Leiden method.",
        epilog="Example: python network_cluster.py in.dist prefix",
    )

    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="cluster file from chdit or vsearch containing cluster information.",
    )

    parser.add_argument(
        "-m",
        "--method",
        metavar="METHOD",
        type=str,
        help="Comparison method containing the necessary information for creating a network.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
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
        "mash": read_in_mash
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
    df = pd.read_csv(file_in, sep="\t", index_col="#query")
    # wide to long
    long_df = df.reset_index().melt(id_vars="#query", var_name="target", value_name="weight")
    # to a network igraph object
    G = ig.Graph.TupleList(long_df.itertuples(index=False), directed=True, weights=True)
    return G


def filter_network(network, threshold):
    """
    Filter the network based on the given score
    """
    filtered_network = network.copy()

    # Get a copy of the edges before removal for iteration
    edges_to_remove = [
        (edge.source, edge.target) for edge in filtered_network.es if edge["weight"] >= threshold or edge["weight"] == 0
    ]

    # Remove edges based on the specified conditions
    filtered_network.delete_edges(edges_to_remove)

    return filtered_network


def cluster_network(network):
    """
    Cluster the network based on the given score
    """
    # Partition the network
    partitions = la.find_partition(
        network, partition_type=la.ModularityVertexPartition, n_iterations=-1, seed=42, weights="weight"
    )

    # extract the names of the vertices
    vertices_names = [[network.vs[index]["name"] for index in cluster] for cluster in partitions]

    return partitions, vertices_names


def to_tsv(vertices_names, prefix):
    """
    Write the network to a tsv file
    """
    # Create a list of lists with the indexed vertices
    indexed_vertices = [[name, idx] for idx, names in enumerate(vertices_names) for name in names]

    # Write the indexed data to a TSV file
    with open(f"{prefix}.tsv", "w") as file:
        for line in indexed_vertices:
            file.write(f"{line[0]}\t{line[1]}\n")


def visualize_network(partitions, network, prefix):
    """
    Visualize the network
    """
    # Set the layout of the network
    layout = network.layout("kk")

    # Plot the network
    ig.plot(partitions, target=f"{prefix}.png", layout=layout, vertex_label=network.vs["name"])


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)

    network = read_in_file(args.file_in, args.method)

    # args.score is ANI, mash calculates distances, so we need to invert the score
    network_filtered = filter_network(network, 1 - args.score)

    clusters, vertices_names = cluster_network(network_filtered)

    to_tsv(vertices_names, args.prefix)

    visualize_network(clusters, network_filtered, args.prefix)


if __name__ == "__main__":
    sys.exit(main())
