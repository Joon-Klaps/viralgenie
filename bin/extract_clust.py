#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import gzip
import json
import logging
import re
import sys
from pathlib import Path

logger = logging.getLogger()


class Cluster:
    """
    A cluster contains the centroid sequence, members of the cluster, cluster_size of centroid.
    """

    def __init__(self, cluster_id, centroid, members):
        self.cluster_id = cluster_id
        self.centroid = centroid
        self.members = members
        if members is not None:
            self.cluster_size = len(members)
        else:
            self.cluster_size = 0

    def set_centroid(self, centroid):
        """
        Set the centroid sequence for the cluster.
        """
        self.centroid = centroid

    def __iter__(self):
        yield "cluster_id", self.row
        yield "centroid", self.centroid
        yield "members", self.members
        yield "cluster_size", self.cluster_size

    def __str__(self):
        return f"Cluster {self.cluster_id} with centroid {self.centroid} and {self.cluster_size} members {self.members}"

    def _save_cluster_members(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_members.txt", "w") as file:
            if self.members:
                for member in self.members:
                    file.write(f"{member}\n")
            else:
                file.write(f"\n")

    def _save_cluster_centroid(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_centroid.txt", "w") as file:
            file.write(f"{self.centroid}\n")

    def _save_cluster_json(self, prefix):
        with open(f"{prefix}_{self.cluster_id}_cluster.json", "w") as file:
            json.dump(self, file, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def _to_line(self, prefix):
        return "\t".join(
            [str(prefix), str(self.cluster_id), str(self.centroid), str(self.cluster_size), ",".join(self.members)]
        )


def parse_clusters_chdit(file_in):
    """
    Extract sequence names from cdhit's cluster files.
    """
    clusters = []
    with open(file_in, "r") as file:
        lines = file.readlines()

    current_cluster_id = None
    current_members = []
    current_centroid = None

    for line in lines:
        if line.startswith(">Cluster"):
            # New cluster detected, add previous Cluster object to the list
            if current_cluster_id is not None:
                cluster = Cluster(current_cluster_id, current_centroid, current_members)
                clusters.append(cluster)

            # Extract the cluster cluster_id from the line and reset the members and centroid
            current_cluster_id = line.strip().split()[1]
            current_members = []
            current_centroid = None
        else:
            # Extract the name from the line
            parts = line.strip().split(">")
            member_name = parts[1].split("...")[0]

            # Check if the line indicates the centroid member
            if parts[1].endswith("*"):
                current_centroid = member_name
            else:
                # Add the member to the list of members
                current_members.append(member_name)

    # Create a Cluster object for the last cluster
    if current_cluster_id is not None:
        cluster = Cluster(current_cluster_id, current_centroid, current_members)
        clusters.append(cluster)

    return clusters


def parse_clusters_vsearch(file_in):
    """
    Extract sequence names from vsearch gzipped cluster files.
    """
    clusters = {}  # Dictionary to store clusters {cluster_id: Cluster}

    with gzip.open(file_in, "rt") as file:
        for line in file:
            line = line.strip()

            if line.startswith("C\t"):
                continue  # Skip the centroid line

            parts = line.split("\t")
            # parts = ['member_type', 'cluster_id', 'length', 'ANI', '...', '...', '...', '...', 'member_name', 'centroid_name']
            cluster_id = parts[1]
            member_name = parts[-2]

            # Create a new cluster object if the cluster cluster_id is not present in the dictionary
            if cluster_id not in clusters.keys():
                clusters[cluster_id] = Cluster(cluster_id, None, [])

            # Set the centroid of the corresponding cluster
            if line.startswith("S\t"):
                clusters[cluster_id].set_centroid(member_name)

            # Append the member to the corresponding cluster
            elif line.startswith("H\t"):
                clusters[cluster_id].members.append(member_name)

    # Convert the dictionary values to a list of clusters and return
    return list(clusters.values())


def print_clusters(clusters, prefix):
    for cluster in clusters:
        cluster._save_cluster_members(prefix)
        cluster._save_cluster_centroid(prefix)
        cluster._save_cluster_json(prefix)

    write_clusters_to_tsv(clusters, prefix)
    write_clusters_summary(clusters, prefix)


def write_clusters_to_tsv(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    with open(f"{prefix}.clusters.tsv", "w") as file:
        file.write("\t".join(["sample", "cluster_id", "centroid", "size", "members"]))
        file.write("\n")
        for cluster in clusters:
            file.write(cluster._to_line(prefix))
            file.write("\n")


def write_clusters_summary(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    avg_size = sum([cluster.cluster_size for cluster in clusters]) / len(clusters)
    n_clusters = len(clusters)
    n_singletons = len([cluster for cluster in clusters if cluster.cluster_size == 0])

    with open(f"{prefix}.summary_mqc.tsv", "w") as file:
        file.write(
            "\n".join(
                [
                    "# id: 'clusters_summary'",
                    "# section_name: 'Clusters summary'",
                    "# description: 'Summary of clusters, displaying the number of clusters, average cluster size and number of singletons (clusters with no members, only a centroid).'",
                    "# format: 'tsv'",
                    "# plot_type: 'table'",
                    "# pconfig:",
                    "#    id: 'clusters_summary'",
                    "#    table_title: 'Clusters summary'",
                ]
            )
        )
        file.write("\n")
        file.write("\t".join(["sample", "n_clusters", "avg_size", "n_singletons"]))
        file.write("\n")
        file.write("\t".join([str(prefix), str(n_clusters), str(avg_size), str(n_singletons)]))
        file.write("\n")


def filter_clusters(clusters, pattern):
    """
    Filter clusters on members given regex pattern.
    """
    filtered_clusters = []
    regex = re.compile(pattern)

    for cluster in clusters:
        if cluster.members:
            matching_members = [member for member in cluster.members if regex.search(member)]
            if matching_members or regex.search(cluster.centroid):
                filtered_clusters.append(Cluster(cluster.cluster_id, cluster.centroid, matching_members))
        elif regex.search(cluster.centroid):
            filtered_clusters.append(cluster)
    return filtered_clusters


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to extract sequence names from cdhit's cluster files.",
        epilog="Example: python extract_cluster.py [cdhitest|vsearch] in.clstr prefix",
    )
    parser.add_argument(
        "option",
        metavar="OPTION",
        type=str,
        help=".clstr file from cdhitestcontaining cluster information.",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="cluster file from chdit or vsearch containing cluster information.",
    )
    parser.add_argument(
        "file_out_prefix",
        metavar="FILE_OUT_PREFIX",
        type=str,
        help="Output file prefix",
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


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)

    if args.option == "cdhitest":
        cluster_list = parse_clusters_chdit(args.file_in)
    elif args.option == "vsearch":
        cluster_list = parse_clusters_vsearch(args.file_in)
    else:
        logger.error(f"Option {args.option} is not supported!")
        sys.exit(2)
    filtered_clusters = filter_clusters(cluster_list, args.pattern)

    print_clusters(filtered_clusters, args.file_out_prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
