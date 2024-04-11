#!/usr/bin/env python

"""Cluster based on a the created network generated from tools"""

import argparse
import logging
import sys
from pathlib import Path

import igraph as ig
import leidenalg as la
import pandas as pd

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
        "-a",
        "--cluster-algorithm",
        metavar="CLUSTER-ALGORITHM",
        type=str,
        default="connected_components",
        help="Algorithm to use for clustering.",
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
        "-c",
        "--chunksize",
        metavar="CHUNKSIZE",
        help="The chunksize to read in the dataframe",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default INFO).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def read_in_file(args):
    """
    Read in the file and return a networkx graph object
    """
    method = args.method
    method_dict = {
        "mash": read_in_mash
        # Add more methods here if needed
    }

    # Check if the method exists in the dictionary
    if method in method_dict:
        logger.info("Choosing method %s", method_dict[method])
        return method_dict[method](args)
    else:
        raise ValueError(f"Method '{method}' not found.")


def read_in_mash(args):
    """
    Read in the file and return a networkx graph object
    """
    INPUT= args.file_in
    CHUNKSIZE= args.chunksize
    THRESHOLD= 1 - args.score # args.score is ANI, mash calculates distances, so we need to invert the score

    output = pd.DataFrame()
    logger.info("Read in file %s", INPUT)
    # see issue 105
    with pd.read_csv(INPUT, sep="\t", encoding="utf-8", index_col="#query", chunksize=CHUNKSIZE) as reader:
        for chunk in reader:
            # wide to long
            long_df = chunk.reset_index().melt(id_vars="#query", var_name="target", value_name="weight")

            # Select only lower triangle
            lower_triangle = long_df[long_df['#query'] >= long_df['target']]

            # Remove what we don't need
            filtered_df = lower_triangle[(lower_triangle['weight'] < THRESHOLD) & (lower_triangle['weight'] != 0)]

            output = pd.concat([output,filtered_df], ignore_index=True)

            logger.info("Processed chunk, current dataframe size: %d", len(output.index))

    # Create Igraph object
    G = ig.Graph.TupleList(output.itertuples(index=False), directed=False, weights=True)
    logger.info("Created the network graph with %d nodes", len(G.vs))

    return G

def cluster_network(network, method):
    """
    Cluster the network based on the given score
    """

    if method == "leiden":
        # Partition the network
        partitions = la.find_partition(
            network, partition_type=la.ModularityVertexPartition, n_iterations=-1, seed=42, weights="weight"
        )

    elif method == "connected_components":
        partitions = network.components(mode="weak")
    else:
        raise ValueError(f"Method '{method}' not found.")
    logger.info("Partitioned the network using %s", method )

    # extract the names of the vertices
    vertices_names = [[network.vs[index]["name"] for index in cluster] for cluster in partitions]
    logger.info("Extracted members of network groups")

    return partitions, vertices_names


def to_tsv(vertices_names, prefix):
    """
    Write the network to a tsv file
    """
    # Create a list of lists with the indexed vertices
    indexed_vertices = [[name, idx] for idx, names in enumerate(vertices_names) for name in names]

    logger.info("Writing network to file")

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
    logger.info("Determined layout of network")

    # Plot the network
    ig.plot(partitions, target=f"{prefix}.png", layout=layout, vertex_label=network.vs["name"])
    logger.info("Visualised the network")



def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    logger.info("Start clustering")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)

    network = read_in_file(args)

    clusters, vertices_names = cluster_network(network, args.cluster_algorithm)

    to_tsv(vertices_names, args.prefix)

    visualize_network(clusters, network, args.prefix)

    logger.info("All done!")


if __name__ == "__main__":
    sys.exit(main())
