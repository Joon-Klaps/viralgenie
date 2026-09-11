#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import gzip
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from utils.file_tools import index_fasta

logger = logging.getLogger()

# Global variables
PATTERN = "^(TRINITY)|(NODE)|(k\d+)|(scaffold\d+)"

# Pre-compiled regex patterns for assembler detection (used in sum_dict_values_by_assembler)
ASSEMBLER_PATTERNS = {
    "spades": re.compile(r"^(spades)|(NODE)"),
    "megahit": re.compile(r"^(megahit)|(k\d+)"),
    "trinity": re.compile(r"^(trinity)|(TRINITY)"),
}

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
        self.cumulative_read_depth = {}

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
            json.dump(
                self,
                file,
                default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__,
                sort_keys=True,
                indent=4,
            )

    def save_centroid_fasta(self, seq_index, prefix):
        """
        Extract the centroid sequence using memory-efficient processing.

        Args:
            seq_index: Pre-built SeqIO index
            prefix: Output file prefix
        """
        centroid_id = self.centroid.split(" ")[0]
        with open(f"{prefix}_{self.cluster_id}_centroid.fa", "w") as file:
            if centroid_id in seq_index:
                SeqIO.write(seq_index[centroid_id], file, "fasta")
            else:
                logger.warning("Centroid %s not found in sequence index", centroid_id)

    def save_members_fasta(self, seq_index, prefix):
        """
        Extract member sequences using pre-built index.

        Args:
            seq_index: Pre-built SeqIO index
            prefix: Output file prefix
        """
        with open(f"{prefix}_{self.cluster_id}_members.fa", "w") as file:
            if not self.members:
                file.write("\n")
                return

            for member in self.members:
                member_id = member.split(" ")[0]
                if member_id in seq_index:
                    SeqIO.write(seq_index[member_id], file, "fasta")
                else:
                    logger.warning("Member %s not found in sequence index", member_id)

    def _to_line(self, prefix):
        rounded_depth = np.round(list(self.cumulative_read_depth.values()), 2).tolist()
        return "\t".join(
            [
                str(prefix),
                str(self.cluster_id),
                str(self.taxid),
                str(self.centroid),
                str(self.cluster_size),
                "\t".join(map(str, rounded_depth)),
                ",".join(self.members),
            ]
        )

    def determine_cumulative_read_depth(self, coverages: Dict) -> Dict:
        """
        Determine the cumulative read depth for each member of the cluster.
        """
        self.cumulative_read_depth = sum_dict_values_by_assembler(self.centroid, self.members, coverages)


def parse_clusters_chdit(file_in: Path, **kwargs) -> List["Cluster"]:
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

    return clusters

def parse_clusters_mmseqs(file_in: Path, **kwargs) -> List["Cluster"]:
    """
    Extract sequence names from mmseqs createtsv output.

    Format:
    Centroid\tMember1\n
    Centroid\tMember2\n
    ...

    """
    taxid = get_taxid(file_in)

    # Read TSV in one shot with pandas - this is highly optimized C code
    df = pd.read_csv(file_in, sep="\t", header=None, names=["centroid", "member"], dtype=str)

    # Group by centroid and aggregate members into lists - vectorized operation
    grouped = df.groupby("centroid", sort=False)["member"].apply(list).reset_index()

    # Assign cluster IDs based on order of first occurrence
    clusters = []
    for idx, row in enumerate(grouped.itertuples(index=False)):
        centroid_name = row.centroid
        members = row.member
        # Members list includes the centroid itself in mmseqs output, filter it out
        members_filtered = [m for m in members if m != centroid_name]
        cluster = Cluster(f"cl{idx}", centroid_name, members_filtered, taxid=taxid)
        clusters.append(cluster)

    return clusters

def parse_clusters_vsearch(file_in: Path, **kwargs) -> List["Cluster"]:
    """
    Extract sequence names from vsearch gzipped cluster files.

    Format:
    record_type\tcluster_id\tlength\tidentity\tstrand\tunused\tunused\tcigar\tquery\ttarget\n
    where record_type is one of S (seed), H (hit), or C (cluster
    for cluster summary).

    """
    taxid = get_taxid(file_in)

    # Read gzipped file with pandas
    # UC format: record_type, cluster_id, length, identity, strand, unused, unused, cigar, query, target
    df = pd.read_csv(
        file_in,
        sep="\t",
        header=None,
        names=["type", "cluster_id", "length", "identity", "strand", "col5", "col6", "cigar", "query", "target"],
        dtype={"type": str, "cluster_id": str, "query": str, "target": str},
        compression="gzip",
        usecols=["type", "cluster_id", "query"],  # Only read needed columns
    )

    # Filter out C (cluster) records - we only need S (seed) and H (hit)
    df = df[df["type"] != "C"].copy()

    # Extract centroids (S records) and members (H records)
    centroids_df = df[df["type"] == "S"][["cluster_id", "query"]].rename(columns={"query": "centroid"})
    members_df = df[df["type"] == "H"][["cluster_id", "query"]].rename(columns={"query": "member"})

    # Group members by cluster_id
    members_grouped = members_df.groupby("cluster_id", sort=False)["member"].apply(list).reset_index()

    # Merge centroids with grouped members
    result = centroids_df.merge(members_grouped, on="cluster_id", how="left")

    # Build Cluster objects
    clusters = []
    for row in result.itertuples(index=False):
        cluster_id = f"cl{row.cluster_id}"
        members = row.member if isinstance(row.member, list) else []
        cluster = Cluster(cluster_id, row.centroid, members, taxid=taxid)
        clusters.append(cluster)

    return clusters

def parse_clusters_clusty(file_in: Path, skip_header: bool = True, **kwargs) -> List["Cluster"]:
    """
    Extract sequence names from clustly cluster files.
    Format:
    sequence_name\tcluster_representative(\tother_info)*

    """
    pattern_regex = re.compile(kwargs.get("pattern", PATTERN))
    taxid = get_taxid(file_in)

    # Read with pandas - handles header skipping and parsing efficiently
    df = pd.read_csv(
        file_in,
        sep="\t",
        header=0 if skip_header else None,
        usecols=[0, 1],  # Only need first two columns
        dtype=str,
    )
    df.columns = ["sequence", "representative"]

    # Extract first word from sequence names (split by whitespace)
    df["sequence"] = df["sequence"].str.split().str[0]

    # Group by representative and aggregate sequences into lists
    grouped = df.groupby("representative", sort=False)["sequence"].apply(list).reset_index()

    # Build clusters
    clusters = []
    for idx, row in enumerate(grouped.itertuples(index=False)):
        key = row.representative
        values = row.sequence

        # Determine centroid
        centroid = key if key in values else get_first_not_match(pattern_regex, values)
        members = [v for v in values if v != centroid]

        cluster = Cluster(f"cl{idx}", centroid, members, taxid=taxid)
        clusters.append(cluster)

    return clusters


def get_taxid(file_in: Path) -> Optional[str]:
    """Extract taxid from file name."""
    pattern = r"_taxid(\d+)_"
    match = re.search(pattern, file_in.name)
    if match:
        logger.debug("Reading %s with taxon id %s...", file_in, match.group(1))
        return match.group(1)
    else:
        logger.debug("No taxon id found for %s", file_in)
        return None

def get_first_not_match(regex_pattern: re.Pattern, data_list: List[str]) -> str:
    """Return the first element that does NOT match the regex_pattern, else return the first element."""
    for item in data_list:
        if not regex_pattern.search(item):
            return item
    return data_list[0]

def write_clusters(clusters: List['Cluster'], sequences: Path, prefix: str, length_clusters: int) -> None:
    """
    Write the clusters to a fasta, json, tsv file.
    """
    logger.info("Building sequence index from %s (this may take a moment for large files)...", sequences)
    seq_index = index_fasta(sequences)
    logger.info("Sequence index built with %d sequences.", len(seq_index))

    for cluster in clusters:
        cluster.save_cluster_members(prefix)
        cluster.save_cluster_centroid(prefix)
        cluster.save_centroid_fasta(seq_index, prefix)
        cluster.save_members_fasta(seq_index, prefix)
        cluster.save_cluster_json(prefix)

    if hasattr(seq_index, 'close'):
        seq_index.close()
    logger.info("Finished writing %d clusters.", len(clusters))
    write_clusters_summary(clusters, prefix, length_clusters)

def write_clusters_to_tsv(clusters, prefix):
    """
    Write the clusters to a json file.
    """
    with open(f"{prefix}.clusters.tsv", "w") as file:
        assemblers = [f"cumulative read depth - {assembler} [%]" for assembler in clusters[0].cumulative_read_depth.keys()]
        file.write("\t".join(["sample", "cluster", "taxon-id", "centroid", "number of members"] + assemblers + ["members"]))
        file.write("\n")
        for cluster in clusters:
            file.write(cluster._to_line(prefix))
            file.write("\n")

def read_coverages(coverages):
    """
    Read the coverages from each idxstats file and compute the percentage of each coverage.
    Return a list of dictionaries, one for each file.

    Format:
    name\tlength\tmapped

    """
    all_coverages = []

    for coverage_file in coverages:
        # Read idxstats file with pandas - typically has 4 columns: name, length, mapped, unmapped
        df = pd.read_csv(
            coverage_file,
            sep="\t",
            header=None,
            usecols=[0, 2],  # Only need name and coverage columns
            names=["name", "coverage"],
            dtype={"name": str, "coverage": np.int64},
        )

        # Aggregate by name (in case of duplicates) and compute percentages vectorized
        coverage_sum = df.groupby("name", sort=False)["coverage"].sum()
        total_coverage = coverage_sum.sum()

        if total_coverage > 0:
            coverage_pct = (coverage_sum / total_coverage) * 100
        else:
            coverage_pct = coverage_sum

        all_coverages.append(coverage_pct.to_dict())

    return all_coverages

def update_cluster_ids(clusters: List["Cluster"]) -> List["Cluster"]:
    # Sort clusters by centroid name before assigning IDs
    sorted_clusters = sorted(clusters, key=lambda c: c.centroid or "")
    updated_clusters = []
    for idx, cluster in enumerate(sorted_clusters):
        cluster.set_cluster_id(f"cl{idx}")
        updated_clusters.append(cluster)
    return updated_clusters

def _get_assembler(key: str) -> Optional[str]:
    """
    Get the assembler type for a given key using pre-compiled patterns.

    Returns:
        Assembler name or None if no pattern matches
    """
    for assembler, pattern in ASSEMBLER_PATTERNS.items():
        if pattern.match(key):
            return assembler
    return None


def sum_dict_values_by_assembler(centroid: str, members: List, coverages: List) -> Dict:
    """
    Sum the values of the dictionaries in the list for the given key.

    Parameters:
    centroid (str): The key to sum the values for.
    members (list): The list of keys to sum the values for.
    coverages (list): The list of dictionaries to sum the values from.

    Returns:
    dict: A dictionary with the summed values for each assembler

    Uses numpy for vectorized summation across coverage dictionaries.
    """
    # Initialize result dictionary using pre-compiled patterns
    result = {assembler: 0.0 for assembler in ASSEMBLER_PATTERNS}

    # Combine all keys to process
    all_keys = [centroid] + members

    # Pre-classify all keys by assembler to avoid repeated regex matching
    assembler_keys = {assembler: [] for assembler in ASSEMBLER_PATTERNS}
    for key in all_keys:
        assembler = _get_assembler(key)
        if assembler is not None:
            assembler_keys[assembler].append(key)

    # For each assembler, sum coverages across all its keys using numpy
    for assembler, keys in assembler_keys.items():
        if keys:
            # Build coverage matrix for vectorized sum
            coverage_values = np.array(
                [[d.get(key, 0.0) for d in coverages] for key in keys],
                dtype=np.float64
            )
            result[assembler] = np.sum(coverage_values)

    return result

def write_clusters_summary(clusters: List["Cluster"], prefix:str, length_clusters:int =None)-> None:
    """
    Write the clusters to a json file.
    """
    avg_size = sum([cluster.cluster_size for cluster in clusters]) / len(clusters)
    n_clusters = len(clusters)
    n_singletons = len([cluster for cluster in clusters if cluster.cluster_size == 0])

    with open(f"{prefix}.summary_mqc.tsv", "w") as file:
        if length_clusters:
            file.write("\t".join(["Sample name", "# Clusters", "# Removed clusters", "Average cluster size", "Number of singletons"]))
            file.write("\n")
            file.write("\t".join([str(prefix), str(n_clusters), str(length_clusters - n_clusters), str(avg_size), str(n_singletons)]))
            file.write("\n")
        else:
            file.write("\t".join(["Sample name", "# Clusters", "Average cluster size", "Number of singletons"]))
            file.write("\n")
            file.write("\t".join([str(prefix), str(n_clusters), str(avg_size), str(n_singletons)]))
            file.write("\n")

def write_clusters_mqc(clusters: List["Cluster"], prefix:str)-> None:
    """
    Write the clusters to a json file.
    """
    with open(f"{prefix}.clusters_mqc.json", "w") as file:
        json.dump(
            clusters,
            file,
            default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o.__dict__,
            sort_keys=True,
            indent=4,
        )

def filter_members(clusters: List["Cluster"], pattern: str = PATTERN) -> List["Cluster"]:
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

def filter_clusters_by_coverage(
    clusters: List["Cluster"], coverages: Dict, threshold: float, keep_n_clusters: int
) -> Tuple[List["Cluster"], List["Cluster"]]:
    """
    Filter clusters on coverage, only keep clusters with a coverage above the threshold. If no clusters are kept, return top 5.
    """
    filtered_clusters = []
    for cluster in clusters:
        cluster.determine_cumulative_read_depth(coverages)
        logger.debug("Cluster %s has cumulative read depth %s",
            cluster.cluster_id, cluster.cumulative_read_depth)
        cum_coverages = np.array(list(cluster.cumulative_read_depth.values()))
        if any(cum_coverages >= threshold):
            filtered_clusters.append(cluster)

    if filtered_clusters:
        return clusters, filtered_clusters

    # Fix: Sum only the values of the cumulative_read_depth dictionary
    sorted_clusters = sorted(clusters, key=lambda x: sum(x.cumulative_read_depth.values()), reverse=True)
    return sorted_clusters, sorted_clusters[:keep_n_clusters]

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
        default=5,
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
        default=PATTERN,
    )

    parser.add_argument(
        "-n",
        "--keep-clusters",
        metavar="KEEP_CLUSTERS",
        type=int,
        help="Define the number of clusters to keep based if no clusters have at least the threshold value's read depth.",
        default=5,
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
    CLUSTER_PARSERS = {
        "mash":            { "func": parse_clusters_clusty  },
        "vrhyme":          { "func": parse_clusters_clusty  },
        "vsearch":         { "func": parse_clusters_vsearch },
        "cdhitest":        { "func": parse_clusters_chdit   },
        "mmseqs-cluster":  { "func": parse_clusters_mmseqs  },
        "mmseqs-linclust": { "func": parse_clusters_mmseqs  }
    }

    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.seq.is_file():
        logger.error("The given input file %s was not found!", args.seq)
        sys.exit(2)

    cluster_list = []

    parser = CLUSTER_PARSERS.get(args.method)
    if not parser:
        logger.error("Option %s is not supported!", args.method)
        sys.exit(2)

    for cluster_file in args.clusters:
        logger.debug("Reading %s...", cluster_file)
        if not cluster_file.is_file():
            logger.error("The given input file %s was not found!", cluster_file)
            sys.exit(2)
        # Call parser function with consistent interface
        cluster_list += parser["func"](cluster_file, pattern=args.pattern, **parser)

    logger.info("Found %s clusters.", len(cluster_list))

    # Remove clusters with no members and external reference
    filtered_clusters = filter_members(cluster_list, args.pattern)
    logger.info("Filtered clusters by members, %s were removed.",
        len(cluster_list) - len(filtered_clusters))

    # Set external reference, used to know if it needs to collapse or called consensus normally
    logger.info("Setting external reference for clusters.")
    for cluster in filtered_clusters:
        cluster.set_external_reference(args.pattern)

    # redefine cluster ids
    clusters_renamed = update_cluster_ids(filtered_clusters)
    logger.info("Renamed cluster ids.")

    clusters = clusters_renamed.copy()
    # Filter clusters by coverage
    if args.coverages:
        coverages = read_coverages(args.coverages)
        clusters, clusters_renamed = filter_clusters_by_coverage(clusters_renamed, coverages, args.perc_reads_contig, args.keep_clusters)
        logger.info("Filtered clusters by coverage, %s were removed.",
            len(filtered_clusters) - len(clusters_renamed))

    assert len(clusters_renamed) != 0, "No clusters left after filtering."

    # Write the clusters to files
    logger.info("Writing results to files.")
    write_clusters_to_tsv(clusters, args.prefix)
    write_clusters(clusters_renamed, args.seq, args.prefix, len(clusters))

    return 0

if __name__ == "__main__":
    sys.exit(main())
