#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import re
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class Cluster:
    """
    A Cdhit cluster contains the reference sequence, members of the cluster, size of reference.
    """

    def __init__(self, id, reference, members):
        self.id = id
        self.reference = reference
        self.members = members

    def __iter__(self):
        yield "id", self.row
        yield "reference", self.reference
        yield "members", self.members

    def __str__(self):
        return f"Cluster {self.id} with reference {self.reference} and members {self.members}"

    def _save_cluster_members(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.id}_members.txt", "w") as file:
            if self.members:
                for member in self.members:
                    file.write(f"{member}\n")
            else:
                file.write(f"\n")

    def _save_cluster_reference(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.id}_reference.txt", "w") as file:
            file.write(f"{self.reference}\n")


def parse_clusters(file_in):
    """
    Extract sequence names from cdhit's cluster files.
    """
    clusters = []
    with open(file_in, "r") as file:
        lines = file.readlines()

    current_cluster_id = None
    current_members = []
    current_reference = None

    for line in lines:
        if line.startswith(">Cluster"):
            # New cluster detected, add previous Cluster object to the list
            if current_cluster_id is not None:
                cluster = Cluster(current_cluster_id, current_reference, current_members)
                clusters.append(cluster)

            # Extract the cluster ID from the line and reset the members and reference
            current_cluster_id = line.strip().split()[1]
            current_members = []
            current_reference = None
        else:
            # Extract the name from the line
            parts = line.strip().split(">")
            member_name = parts[1].split("...")[0]

            # Check if the line indicates the reference member
            if parts[1].endswith("*"):
                current_reference = member_name
            else:
                # Add the member to the list of members
                current_members.append(member_name)

    # Create a Cluster object for the last cluster
    if current_cluster_id is not None:
        cluster = Cluster(current_cluster_id, current_reference, current_members)
        clusters.append(cluster)

    return clusters


def filter_clusters(clusters, pattern):
    """
    Filter clusters on members given regex pattern.
    """
    filtered_clusters = []
    regex = re.compile(pattern)

    for cluster in clusters:
        if cluster.members:
            all_members = [cluster.reference] + cluster.members
        else:
            all_members = [cluster.reference]

        # Check if any members match the pattern
        matching_members = [member for member in all_members if regex.search(member)]

        # If there are matching members, keep the cluster
        if matching_members:
            filtered_clusters.append(Cluster(cluster.id, cluster.reference, matching_members))

    return filtered_clusters


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to extract sequence names from cdhit's cluster files.",
        epilog="Example: python extract_cdhit.py chdit.clstr prefix",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help=".clstr file from cdhit containing cluster information.",
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
        help="Regex pattern to filter clusters by reference sequence name.",
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

    cluster_list = parse_clusters(args.file_in)
    filtered_clusters = filter_clusters(cluster_list, args.pattern)

    for cluster in filtered_clusters:
        cluster._save_cluster_members(args.file_out_prefix)
        cluster._save_cluster_reference(args.file_out_prefix)


if __name__ == "__main__":
    sys.exit(main())
