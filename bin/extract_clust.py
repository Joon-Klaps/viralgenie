#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import gzip
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger()


class Cluster:
    """
    A cluster contains the centroid sequence, members of the cluster, cluster_size of centroid.
    """

    def __init__(self, cluster_id, centroid, members):
        # Having only a number get's removed by multiqc which causes merging errors downstream
        self.cluster_id = cluster_id
        self.centroid = centroid
        self.members = members
        self.external_reference = None
        if members is not None:
            self.cluster_size = len(members)
        else:
            self.cluster_size = 0

    def _set_centroid(self, centroid):
        """
        Set the centroid sequence for the cluster.
        """
        self.centroid = centroid

    def _set_cluster_id(self, id):
        """
        Set the centroid sequence for the cluster.
        """
        self.cluster_id = id

    def set_external_reference(self, pattern):
        """
        Set the external reference for the cluster.
        """
        regex = re.compile(pattern)
        self.external_reference = not bool(regex.search(self.centroid))

    def __iter__(self):
        yield "cluster_id", self.row
        yield "centroid", self.centroid
        yield "members", self.members
        yield "cluster_size", self.cluster_size

    def __str__(self):
        return f"Cluster {self.cluster_id} with centroid {self.centroid}, external {self.external_reference} and {self.cluster_size} members {self.members}"

    def _save_cluster_members(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_members.txt", "w") as file:
            if self.members:
                for member in self.members:
                    file.write(f"{member}\n")
            else:
                file.write(f"")

    def _save_cluster_centroid(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_centroid.txt", "w") as file:
            file.write(f"{self.centroid}\n")

    def _save_cluster_json(self, prefix):
        with open(f"{prefix}_{self.cluster_id}_cluster.json", "w") as file:
            json.dump(self, file, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def _save_centroid_fasta(self, sequences, prefix):
        """
        Extract the sequences from the input file based on the groups.
        """
        sequence_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
        with open(f"{prefix}_{self.cluster_id}_centroid.fa", "w") as file:
            centroid_id = self.centroid.split(" ")[0]
            if centroid_id in sequence_dict:
                SeqIO.write(sequence_dict[centroid_id], file, "fasta")

    def _save_members_fasta(self, sequences, prefix):
        """
        Extract the sequences from the input file based on the groups.
        """
        sequence_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
        with open(f"{prefix}_{self.cluster_id}_members.fa", "w") as file:
            if self.members:
                for member in self.members:
                    seq_name = member.split(" ")[0]
                    if seq_name in sequence_dict:
                        SeqIO.write(sequence_dict[seq_name], file, "fasta")
            else:
                file.write(f"\n")

    def _to_line(self, prefix):
        return "\t".join(
            [str(prefix), str(self.cluster_id), str(self.centroid), str(self.cluster_size), ",".join(self.members)]
        )


def parse_clusters_chdit(file_in):
    """
    Extract sequence names from cdhit's cluster files.
    """
    with open(file_in, "r") as file:
        lines = file.readlines()

    clusters = []
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
            current_cluster_id = f"cl{line.strip().split()[1]}"
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

    return clusters.copy()


def parse_clusters_mmseqs(file_in):
    """
    Extract sequence names from mmseqs createtsv output.
    """

    # Dictionary to store clusters {cluster_id: Cluster}
    clusters = {}

    with open(file_in, "rt") as file:
        for line in file:
            centroid_name, member_name = line.strip().split("\t")
            cluster_id = f"cl{len(clusters)}"  # Generate unique cluster ID

            # Check if centroid already exists in clusters
            existing_cluster = next((c for c in clusters.values() if c.centroid == centroid_name), None)

            if existing_cluster:
                existing_cluster.members.append(member_name)
                existing_cluster.cluster_size = len(existing_cluster.members)
            else:
                new_cluster = Cluster(cluster_id, centroid_name, [])
                clusters[cluster_id] = new_cluster
    return list(clusters.values())


def parse_clusters_vsearch(file_in):
    """
    Extract sequence names from vsearch gzipped cluster files.
    """
    # Dictionary to store clusters {cluster_id: Cluster}
    clusters = {}

    with gzip.open(file_in, "rt") as file:
        for line in file:
            line = line.strip()

            if line.startswith("C\t"):
                continue  # Skip the centroid line

            parts = line.split("\t")
            # parts = ['member_type', 'cluster_id', 'length', 'ANI', '...', '...', '...', '...', 'member_name', 'centroid_name']
            cluster_id = f"cl{parts[1]}"
            member_name = parts[-2]

            # Create a new cluster object if the cluster cluster_id is not present in the dictionary
            if cluster_id not in clusters.keys():
                clusters[cluster_id] = Cluster(cluster_id, None, [])

            # Set the centroid of the corresponding cluster
            if line.startswith("S\t"):
                clusters[cluster_id]._set_centroid(member_name)

            # Append the member to the corresponding cluster
            elif line.startswith("H\t"):
                clusters[cluster_id].members.append(member_name)

    # Convert the dictionary values to a list of clusters and return
    return list(clusters.values())


def parse_clusters_vrhyme(file_in, pattern, skip_header=True):
    """
    Extract sequence names from vrhyme gzipped cluster files using regex.
    input file:
        scaffold	bin
        k39_0 flag=1 multi=111.6808 len=29849	1
        Japan/DP0078/2020	1
        CHN/HN04/2020	1
    """
    clusters = {}  # Dictionary to store clusters {cluster_id: Cluster}
    grouped = defaultdict(list)
    pattern_regex = re.compile(pattern)

    with open(file_in, "r") as file:
        if skip_header:
            next(file)
        for line in file:
            value, key = line.strip().split("\t")
            grouped[key].append(value.split()[0])

        for key, values in grouped.items():
            centroid = get_first_not_match(pattern_regex, values)
            members = [value for value in values if value != centroid]
            cluster_id = f"cl{key}"
            clusters[cluster_id] = Cluster(cluster_id, centroid, members)

    return list(clusters.values())


def get_first_not_match(regex_pattern, data_list):
    """
    Return the first element that matches the regex_pattern else return the first element.
    """
    for item in data_list:
        match = re.search(regex_pattern, item)
        if not match:
            return item
    return data_list[0]


def write_clusters(clusters, sequences, prefix):
    for cluster in clusters:
        cluster._save_cluster_members(prefix)
        cluster._save_cluster_centroid(prefix)
        cluster._save_centroid_fasta(sequences, prefix)
        cluster._save_members_fasta(sequences, prefix)
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


def update_cluster_ids(clusters):
    updated_clusters = []
    for idx, cluster in enumerate(clusters):
        cluster._set_cluster_id(f"cl{idx}")
        updated_clusters.append(cluster)
    return updated_clusters


def write_clusters_summary(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    avg_size = sum([cluster.cluster_size for cluster in clusters]) / len(clusters)
    n_clusters = len(clusters)
    n_singletons = len([cluster for cluster in clusters if cluster.cluster_size == 0])

    with open(f"{prefix}.summary_mqc.tsv", "w") as file:
        file.write("\t".join(["Sample name", "# Clusters", "Average cluster size", "Number of singletons"]))
        file.write("\n")
        file.write("\t".join([str(prefix), str(n_clusters), str(avg_size), str(n_singletons)]))
        file.write("\n")


def filter_clusters(clusters, pattern):
    """
    Filter clusters on members given regex pattern, members cannot contain the pattern.
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
        epilog="Example: python extract_cluster.py [cdhitest|vsearch] --clusters in.clstr1 in.clstr2 ... --seq in.seq prefix",
    )
    parser.add_argument(
        "-m",
        "--method",
        metavar="CLUSTER METHOD",
        type=str,
        choices=("cdhitest", "vsearch", "mmseqs-linclust", "mmseqs-cluster", "mash", "vrhyme"),
        help="Cluster algorithm used to generate cluster files.",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        nargs="+",
        metavar="CLUSTERS_IN",
        type=Path,
        help="cluster file from cluster methods containing cluster information.",
    )
    parser.add_argument(
        "-s",
        "--seq",
        metavar="SEQ_IN",
        type=Path,
        help="cluster file from chdit, vsearch, mmseqs_createtsv containing cluster information.",
    )
    parser.add_argument(
        "-p" "--prefix",
        dest="prefix",
        metavar="PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "-r",
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

    if not args.seq.is_file():
        logger.error(f"The given input file {args.seq} was not found!")
        sys.exit(2)

    cluster_list = []
    for cluster_file in args.clusters:
        if not cluster_file.is_file():
            logger.error(f"The given input file {cluster_file} was not found!")
            sys.exit(2)
        if args.method == "cdhitest":
            cluster_list += parse_clusters_chdit(cluster_file)
        elif args.method == "vsearch":
            cluster_list += parse_clusters_vsearch(cluster_file)
        elif args.method == "mmseqs-linclust" or args.method == "mmseqs-cluster":
            cluster_list += parse_clusters_mmseqs(cluster_file)
        elif args.method == "vrhyme":
            # vrhyme doens't select centroids so we provide pattern to give preference to non matching sequences
            cluster_list += parse_clusters_vrhyme(cluster_file, args.pattern)
        elif args.method == "mash":
            cluster_list += parse_clusters_vrhyme(cluster_file, args.pattern, False)
        else:
            logger.error(f"Option {args.method} is not supported!")
            sys.exit(2)

    # redefine cluster ids
    clusters_renamed = update_cluster_ids(cluster_list)

    # Remove clusters with no members and external reference
    filtered_clusters = filter_clusters(clusters_renamed, args.pattern)

    # Set external reference, used to know if it needs to collapse or called consensus normally
    for cluster in filtered_clusters:
        cluster.set_external_reference(args.pattern)

    write_clusters(filtered_clusters, args.seq, args.prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
