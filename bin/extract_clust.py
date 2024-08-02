#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import gzip
import json
import logging
import re
import sys
import numpy as np
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger()

class Cluster:
    """
    A cluster contains the centroid sequence, members of the cluster, cluster_size of centroid.
    """

    def __init__(self, cluster_id, centroid, members, taxid=None):
        # Having only a number get's removed by multiqc which causes merging errors downstream
        self.cluster_id = cluster_id
        self.centroid = centroid
        self.members = members
        self.external_reference = None
        self.taxid = taxid
        if members is not None:
            self.cluster_size = len(members)
        else:
            self.cluster_size = 0
        self.cumulative_read_depth = []

    def _set_centroid(self, centroid):
        """
        Set the centroid sequence for the cluster.
        """
        self.centroid = centroid

    def set_cluster_id(self, id):
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
        yield "cluster_id", self.cluster_id
        yield "centroid", self.centroid
        yield "members", self.members
        yield "cluster_size", self.cluster_size
        yield "taxid", self.taxid

    def __str__(self):
        return f"Cluster {self.cluster_id}, taxid {self.taxid} with centroid {self.centroid}, external {self.external_reference} and {self.cluster_size} members {self.members}"

    def save_cluster_members(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_members.txt", "w") as file:
            if self.members:
                for member in self.members:
                    file.write(f"{member}\n")
            else:
                file.write(f"")

    def save_cluster_centroid(self, prefix):
        """
        Save the cluster to a file.
        """
        with open(f"{prefix}_{self.cluster_id}_centroid.txt", "w") as file:
            file.write(f"{self.centroid}\n")

    def save_cluster_json(self, prefix):
        with open(f"{prefix}_{self.cluster_id}_cluster.json", "w") as file:
            json.dump(self, file, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__, sort_keys=True, indent=4)

    def save_centroid_fasta(self, sequences, prefix):
        """
        Extract the sequences from the input file based on the groups.
        """
        sequence_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
        with open(f"{prefix}_{self.cluster_id}_centroid.fa", "w") as file:
            centroid_id = self.centroid.split(" ")[0]
            if centroid_id in sequence_dict:
                SeqIO.write(sequence_dict[centroid_id], file, "fasta")

    def save_members_fasta(self, sequences, prefix):
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
        rounded_depth = np.round(self.cumulative_read_depth, 2).tolist()
        return "\t".join(
            [
            str(prefix),
            str(self.taxid),
            str(self.cluster_id),
            str(self.centroid),
            str(self.cluster_size),
            ",".join(map(str,rounded_depth)),
            ",".join(self.members)
            ]
        )

    def determine_cumulative_read_depth(self, coverages):
        """
        Determine the cumulative read depth for each member of the cluster.
        """
        self.cumulative_read_depth = np.sum([
                [d.get(key, 0) for d in coverages] for key in [self.centroid] + self.members
            ], axis=0)

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
    taxid = get_taxid(file_in)

    for line in lines:
        if line.startswith(">Cluster"):
            # New cluster detected, add previous Cluster object to the list
            if current_cluster_id is not None:
                cluster = Cluster(current_cluster_id, current_centroid, current_members, taxid=taxid)
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
        cluster = Cluster(current_cluster_id, current_centroid, current_members, taxid=taxid)
        clusters.append(cluster)

    return clusters.copy()

def parse_clusters_mmseqs(file_in):
    """
    Extract sequence names from mmseqs createtsv output.
    """

    # Dictionary to store clusters {cluster_id: Cluster}
    clusters = {}
    taxid = get_taxid(file_in)

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
                new_cluster = Cluster(cluster_id, centroid_name, [], taxid=taxid)
                clusters[cluster_id] = new_cluster
    return list(clusters.values())

def parse_clusters_vsearch(file_in):
    """
    Extract sequence names from vsearch gzipped cluster files.
    """
    # Dictionary to store clusters {cluster_id: Cluster}
    clusters = {}
    taxid = get_taxid(file_in)

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
                clusters[cluster_id] = Cluster(cluster_id, None, [], taxid=taxid)

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
    taxid = get_taxid(file_in)

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
            clusters[cluster_id] = Cluster(cluster_id, centroid, members, taxid=taxid)

    return list(clusters.values())

def get_taxid(file_in):
    """
    Extract taxid from file name.

    Parameters:
    file_in (str or file-like object): The file name or file object from which to extract the taxid.

    Returns:
    str or None: The extracted taxid if found, None otherwise.
    """

    pattern = r'_taxid(\d+)_'
    match = re.search(pattern, file_in.name)
    if match:
        logger.debug(f"Reading {file_in} with taxon id {match.group(1)}...")
        return match.group(1)
    else :
        logger.debug(f"No taxon id found for {file_in} ")
        return None

def get_first_not_match(regex_pattern, data_list):
    """
    Return the first element that matches the regex_pattern else return the first element.
    """
    for item in data_list:
        match = re.search(regex_pattern, item)
        if not match:
            return item
    return data_list[0]

def write_clusters(clusters, sequences, prefix) -> None:
    """
    Write the clusters to a fasta, json, tsv file.
    """
    for cluster in clusters:
        cluster.save_cluster_members(prefix)
        cluster.save_cluster_centroid(prefix)
        cluster.save_centroid_fasta(sequences, prefix)
        cluster.save_members_fasta(sequences, prefix)
        cluster.save_cluster_json(prefix)

    write_clusters_summary(clusters, prefix)

def write_clusters_to_tsv(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    with open(f"{prefix}.clusters.tsv", "w") as file:
        file.write("\t".join(["sample", "taxon-id", "cluster-id", "centroid", "size","cumulative read depth [%]", "members"]))
        file.write("\n")
        for cluster in clusters:
            file.write(cluster._to_line(prefix))
            file.write("\n")

def read_coverages(coverages):
    """
    Read the coverages from each idxstats file and compute the percentage of each coverage.
    Return a list of dictionaries, one for each file.
    """
    all_coverages = []

    for coverage_file in coverages:
        coverages_dict = {}
        total_coverage = 0

        # First pass to compute the total coverage for this file
        with open(coverage_file, "r") as file:
            for line in file:
                parts = line.strip().split("\t")
                coverage = int(parts[2])
                total_coverage += coverage
                if parts[0] in coverages_dict:
                    coverages_dict[parts[0]] += coverage
                else:
                    coverages_dict[parts[0]] = coverage

        # Compute the percentage for each entry
        for key in coverages_dict:
            coverages_dict[key] = (coverages_dict[key] / total_coverage) * 100

        # Add the coverage dictionary for this file to the list
        all_coverages.append(coverages_dict)

    return all_coverages

def update_cluster_ids(clusters):
    updated_clusters = []
    for idx, cluster in enumerate(clusters):
        cluster.set_cluster_id(f"cl{idx}")
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

def write_clusters_mqc(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    with open(f"{prefix}.clusters_mqc.json", "w") as file:
        json.dump(clusters, file, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__, sort_keys=True, indent=4)

def filter_members(clusters, pattern):
    """
    Filter clusters on members given regex pattern, members cannot contain the pattern.
    """
    filtered_clusters = []
    regex = re.compile(pattern)

    for cluster in clusters:
        if cluster.members:
            matching_members = [member for member in cluster.members if regex.search(member)]
            if matching_members or regex.search(cluster.centroid):
                filtered_clusters.append(Cluster(cluster.cluster_id, cluster.centroid, matching_members, taxid=cluster.taxid))
        elif regex.search(cluster.centroid):
            filtered_clusters.append(cluster)
    return filtered_clusters

def filter_clusters_by_coverage(clusters: list , coverages: dict, threshold: float) -> list:
    """
    Filter clusters on coverage, only keep clusters with a coverage above the threshold.
    """
    filtered_clusters = []
    for cluster in clusters:
        cluster.determine_cumulative_read_depth(coverages)
        logger.debug("Cluster %s has cumulative read depth %s", cluster.cluster_id, cluster.cumulative_read_depth)
        if any(cluster.cumulative_read_depth >= threshold):
            filtered_clusters.append(cluster)
    return clusters,filtered_clusters

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
        "-d",
        "--coverages",
        nargs="+",
        metavar="COVERAGES",
        type=Path,
        help="idxstats file displaying the number of reads mapped to each contig.",
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
        "-t",
        "--perc_reads_contig",
        default= 5,
        metavar="PERC_READS_CONTIG",
        type=float,
        help="Percentage of reads mapped to contig to keep cluster.",
    )

    parser.add_argument(
        "-r",
        "--pattern",
        metavar="PATTERN",
        type=str,
        help="Regex pattern to filter clusters by centroid sequence name.",
        default="^(TRINITY)|(NODE)|(k\d+)|(scaffold\d+)",  # Default pattern matches Trinity, SPADes, MEGAHIT, sspace_basic assembly names
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
        logger.debug(f"Reading {cluster_file}...")
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

    logger.info("Found %d clusters.", len(cluster_list))

    # redefine cluster ids
    clusters_renamed = update_cluster_ids(cluster_list)
    logger.info("Renamed cluster ids.")

    # Remove clusters with no members and external reference
    filtered_clusters = filter_members(clusters_renamed, args.pattern)
    logger.info("Filtered clusters by members, %d were removed.", len(clusters_renamed) - len(filtered_clusters))

    # Set external reference, used to know if it needs to collapse or called consensus normally
    logger.info("Setting external reference for clusters.")
    for cluster in filtered_clusters:
        cluster.set_external_reference(args.pattern)

    # Filter clusters by coverage
    if args.coverages:
        coverages = read_coverages(args.coverages)
        clusters,filtered_clusters = filter_clusters_by_coverage(filtered_clusters, coverages, args.perc_reads_contig)
        logger.info("Filtered clusters by coverage, %d were removed.", len(clusters_renamed) - len(filtered_clusters))

    # Write the clusters to files
    logger.info("Writing results to files.")
    write_clusters_to_tsv(clusters, args.prefix)
    write_clusters(filtered_clusters, args.seq, args.prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
