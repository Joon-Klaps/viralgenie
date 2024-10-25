#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import argparse
import csv
import logging
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple


import pandas as pd
import numpy as np
import multiqc as mqc
from multiqc.plots import bargraph, table
from multiqc.types import Anchor
import yaml


# import plotly.io as pio
# import plotly.express as px

logger = logging.getLogger()

BLAST_COLUMNS = [
    "query",
    "subject",
    "subject title",
    "pident",
    "qlen",
    "slen",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

SUMMARY_COLUMNS = {
    "index": "index",
    "sample name": "sample name",
    "cluster": "cluster",
    "step": "step",
    "species": "species",
    "segment": "segment",
    "definition": "definition",
    "(annotation) species": "(annotation) species",
    "(annotation) segment": "(annotation) segment",
    "(blast) length": "contig length",
    "(annotation) % contig aligned": "% contig aligned",
    "(quast) % N's": "% N's",
    "(bcftools_stats) number of SNPs": "number of SNPs",
    "(multiqc) mosdepth Median read depth": "Median read depth",
}

CONSTRAIN_GENERAL_STATS_COLUMNS = [
    "input_reads",
    "output_reads",
    "number_of_SNPs",
    "CLUSTER: mosdepth.mean_coverage",
    "CLUSTER: mosdepth.min_coverage",
    "CLUSTER: mosdepth.max_coverage",
    "CLUSTER: mosdepth.median_coverage",
    "CLUSTER: mosdepth.1_x_pc",
    "CLUSTER: mosdepth.10_x_pc",
]

GROUPING_COLUMNS = [
    "(annotation) species",
    "(annotation) segment",
    "(annotation) lineage",
]

FILES_OF_INTEREST = {
    "samtools": "multiqc_samtools_stats",
    "umitools": "multiqc_umitools_dedup",
    "multiqc_general_stats": "",
    "picard": "mutliqc_picard_dups",
    "ivar_variants": "",
    "bcftools": "multiqc_bcftools_stats",
}

CLUSTER_HEADERS = {
    "# Filtered clusters": {
        "title": "Filtered # clusters",
        "description": "Number of contig clusters used for further refinement",
        "scale": "Blues",
    },
    "Total # clusters": {
        "title": "Total # clusters",
        "description": "Total number of input contig clusters before filtering",
        "scale": "Blues",
    },
    "# Clusters": {
        "title": "# Clusters",
        "description": "Number of contig clusters used for further refinement ",
        "scale": "Blues",
    },
}

CLUSTER_PCONFIG = {
    "id": "summary_clusters_info",
    "title": "Number of contig clusters",
    "ylab": "# clusters",
    "y_decimals": False,
}


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to combine individual log & summary files which we will pass down to multiqc subsequently.",
        epilog="Example: python custom_multiqc_tables.py --clusters_summary file1,file2,file3,... ",
    )

    def file_choices(choices, fname):
        fname_path = Path(fname)
        if not fname_path.is_file():
            logger.error(f"File '{fname}' does not exist")
            sys.exit(2)
        ext = fname_path.suffix[1:]
        if ext not in choices:
            logger.error(f"File '{fname}' with {ext}, doesn't end with one of {choices}")
            sys.exit(2)
        return fname_path

    parser.add_argument(
        "--multiqc_files",
        metavar="MULTIQC FILES",
        help="Input files for the multiqc module",
        type=Path,
    )

    parser.add_argument(
        "--multiqc_config",
        metavar="MULTIQC CONFIG FILE",
        help="Multiqc config file for report structure & layout",
        type=Path,
    )

    parser.add_argument(
        "--clusters_summary",
        metavar="CLUSTER SUMMARY FILES",
        nargs="+",
        type=Path,
        help="List of cluster summary files from created by the module extract_clust.py",
    )

    parser.add_argument(
        "--annotation_files",
        metavar="Annotation FILES",
        nargs="+",
        help="Blast files for each contig to the annotation database, having the standard outfmt 6",
        type=Path,
    )

    parser.add_argument(
        "--prefix",
        metavar="FILE_OUT_PREFIX",
        type=str,
        help="Output file prefix",
    )

    parser.add_argument(
        "--bed_files",
        metavar="BED_FILES",
        nargs="+",
        help="Bed (coverage) files for each sample",
        type=lambda s: file_choices(("bed", "gz"), s),
    )

    parser.add_argument(
        "--sample_metadata",
        metavar="SAMPLE METADATA",
        help="Sample metadata file containing information on the samples, supported formats: '.csv', '.tsv', '.yaml', '.yml'",
        type=lambda s: file_choices(("csv", "tsv", "yaml", "yml"), s),
    )

    parser.add_argument(
        "--screen_files",
        metavar="MASH SCREEN FILES",
        nargs="+",
        help="Mash screen result of the module SELECT REFERENCE where top hit is outputted in json file format",
        type=lambda s: file_choices(("json"), s),
    )

    parser.add_argument(
        "--mapping_constrains",
        metavar="MAPPING CONSTRAINS",
        help="Mapping constrains file containing information on the sequences that need to be used for mapping against the samples, supported formats: '.csv', '.tsv', '.yaml', '.yml'",
        type=lambda s: file_choices(("csv", "tsv", "yaml", "yml"), s),
    )

    parser.add_argument(
        "--comment_dir",
        metavar="MULTIQC COMMENT DIR",
        help="Directory with the multiqc header files for table annotation that correspond to the different tables that will be created",
        type=Path,
    )

    parser.add_argument(
        "--checkv_files",
        metavar="CHECKV FILES",
        nargs="+",
        help="Checkv summary files for each sample",
        type=Path,
    )
    parser.add_argument(
        "--filter_level",
        metavar="FILTER LEVEL",
        choices=["normal", "strict", "none"],
        default="normal",
        type=str,
        help="Specify how strict the filtering should be, default is normal.",
    )

    parser.add_argument(
        "--clusters_files",
        metavar="CLUSTER FILES",
        nargs="+",
        type=Path,
        help="Cluster files for each sample in table format containing information on the number every cluster of a sample.",
    )

    parser.add_argument(
        "--blast_files",
        metavar="BLAST FILES",
        nargs="+",
        help="Blast files for each contig, having the standard outfmt 6",
        type=Path,
    )

    parser.add_argument(
        "--quast_files",
        metavar="QUAST FILES",
        nargs="+",
        help="Quast summary files for each sample",
        type=Path,
    )

    parser.add_argument(
        "--save_intermediate",
        metavar="SAVE INTERMEDIATE FILES",
        type=bool,
        nargs="?",
        const=True,
        default=False,
    )

    parser.add_argument(
        "--table_headers",
        metavar="TABLE HEADERS",
        help="Yaml file with the table headers for the different tables that will be created",
        type=Path,
    )

    parser.add_argument(
        "--multiqc_dir",
        metavar="MULTIQC DIR",
        help="Multiqc directory where the multiqc files will be used to create the custom tables for multiqc",
        type=Path,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)


def check_file_exists(file, throw_error=True):
    """Check if the given files exist.

    Args:
        file (str): The path to the file to be checked.
        throw_error (bool, optional): Whether to throw an error and exit the program if the file is not found.
            Defaults to True.

    Returns:
        bool: True if the file exists and is not empty, False otherwise.
    """
    if not Path(file).exists():
        if throw_error:
            logger.error("The given input file %s was not found!", file)
            sys.exit(2)
        else:
            logger.warning("The given input file %s was not found!", file)
            return False
    elif not os.stat(file).st_size > 0:
        # logger.warning("The given input file %s is empty, it will not be used!", file)
        return False
    return True


def read_header_file(file_path):
    """
    Read a file and return its content as a list of strings.

    Args:
        file_path (str): The path to the file to be read.

    Returns:
        list: A list of strings representing the content of the file.
    """
    with open(file_path, "r") as file:
        content = file.read().splitlines()
    return content


def concat_table_files(table_files, **kwargs):
    """Concatenate all the cluster summary files into a single dataframe.

    Args:
        table_files (list): List of file paths to be concatenated.
        **kwargs: Additional keyword arguments to be passed to pd.read_csv().

    Returns:
        pd.DataFrame: The concatenated dataframe.
    """
    try:
        valid_dfs = [read_file_to_dataframe(file, **kwargs) for file in table_files if check_file_exists(file)]

        if not valid_dfs:
            logging.warning("Warning concatenating files: %s", table_files)
            logging.warning("No valid files found to concatenate.")
            return pd.DataFrame()

        df = pd.concat(valid_dfs)
        return df

    except ValueError as e:
        logging.warning("Error concatenating files: %s\n%s", table_files, e)
        return pd.DataFrame()


def read_in_quast(table_files):
    """Concatenate all the cluster summary files into a single dataframe.

    Args:
        table_files (list): List of file paths to the cluster summary files.

    Returns:
        pd.DataFrame: A dataframe containing the concatenated data from all the cluster summary files.
    """
    df = pd.DataFrame()
    if table_files:
        for file in table_files:
            if check_file_exists(file, throw_error=False):
                with open(file, "r") as f:
                    d = dict(line.strip().split("\t") for line in f)
                    df = pd.concat([df, pd.DataFrame.from_dict(d, orient="index").T])
    return df


def get_module_selection(table_headers: Path = None) -> Dict:
    """
    Get the files of interest and the columns of interest from the table headers file

    Args:
        table_headers (str): Path to the table headers file

    Returns:
        a dictionary containing the {module:{section:{old_col:new_col}}}, both section and old & new column can be empty
    """
    if not table_headers:
        return FILES_OF_INTEREST

    check_file_exists(table_headers)
    yaml_data = yaml.safe_load(table_headers.read_text())

    return yaml_data


def write_dataframe(df, file, comment):
    """
    Write a pandas DataFrame to a file in TSV format.

    Args:
        df (pd.DataFrame): The DataFrame to be written.
        file (str): The file path to write the DataFrame to.
        comment (list): A list of strings to be written as comments at the beginning of the file.

    Returns:
        None
    """
    if df.empty:
        logger.warning("The DataFrame %s is empty, nothing will be written to the file!", file)
        return
    df_tsv = df.to_csv(sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC)
    with open(file, "w") as f:
        if comment:
            f.write("\n".join(comment))
            f.write("\n")
        f.write(df_tsv)


def filter_and_rename_columns(data: pd.DataFrame, columns: List[Union[str, Dict[str, str]]]) -> pd.DataFrame:
    """
    Filter and rename columns in a DataFrame.

    Args:
        data (pd.DataFrame): Input DataFrame.
        columns (List[Union[str, Dict[str, str]]]): List of column names or rename dictionaries.

    Returns:
        pd.DataFrame: Filtered and renamed DataFrame.
    """
    rename_dict = {}
    keep_columns = []

    if not columns:
        return data

    for col in columns:
        if isinstance(col, str):
            keep_columns.append(col)
        elif isinstance(col, dict):
            old_name, new_name = next(iter(col.items()))
            rename_dict[old_name] = new_name
            keep_columns.append(old_name)

    # Check if all columns exist in the DataFrame
    missing_columns = set(keep_columns) - set(data.columns)
    if missing_columns:
        logger.warning("Columns not found in data: %s\n columns available: %s", missing_columns, data.columns)
        keep_columns = [col for col in keep_columns if col in data.columns]

    return data[keep_columns].rename(columns=rename_dict)


def df_from_tsv(file, **kwargs):
    """
    Read a dataframe from a tsv file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as table:
        df = pd.read_csv(table, sep="\t", **kwargs)
    return df


def df_from_csv(file, **kwargs):
    """
    Read a dataframe from a csv file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as table:
        df = pd.read_csv(table, **kwargs)
    return df


def df_from_yaml(file, **kwargs):
    """
    Read a dataframe from a YAML file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as yaml_file:
        data = yaml.safe_load(yaml_file, **kwargs)
        df = pd.DataFrame(data)
    return df


def df_from_json(file, **kwargs):
    """
    Read a dataframe from a JSON file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as json_file:
        try:
            data = json.load(json_file, **kwargs)
        except json.JSONDecodeError as e:
            logger.warning("Error reading JSON file %s: %s", file, e)
            return pd.DataFrame()

        # Check if 'query' key exists
        if "filename" not in data:
            # Get the filename without path and suffix
            filename = os.path.splitext(os.path.basename(file))[0]
            # Add new key-value pair
            data["filename"] = filename + "_constrain"
        df = pd.DataFrame([data])
    return df


def read_file_to_dataframe(file, **kwargs):
    """
    Read a dataframe from a file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    file_path = Path(file)
    if os.path.getsize(file_path) == 0:
        logger.debug("File is empty %s", file_path)
        return pd.DataFrame()
    if file_path.suffix in [
        ".tsv",
        ".txt",
    ]:  # mqc calls tsv's txts, bed files are gzipped
        df = df_from_tsv(file_path, **kwargs)
    elif file_path.suffix == ".csv":
        df = df_from_csv(file_path, **kwargs)
    elif file_path.suffix in [".yaml", ".yml"]:
        df = df_from_yaml(file_path, **kwargs)
    elif file_path.suffix in [".json"]:
        df = df_from_json(file_path, **kwargs)
    else:
        logger.error(
            "The file format %s is not supported of file %s!",
            file_path.suffix,
            file_path,
        )
        sys.exit(2)
    return df


def join_dataframe(base_df: pd.DataFrame, dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Join multiple DataFrames or a DataFrame and a Series together.

    Args:
        base_df (pd.DataFrame or pd.Series): The base DataFrame or Series to start with.
        dfs (list): A list of DataFrames to be joined with the base DataFrame.

    Returns:
        pd.DataFrame: The joined DataFrame.
    """
    joined_df = base_df.copy()

    for df in dfs:
        if not isinstance(df, pd.DataFrame):
            logger.error("Unable to process the df: %s of class %s", df, type(df))
        if df.empty:
            continue
        if isinstance(joined_df, pd.Series):
            joined_df = joined_df.to_frame()
        joined_df = pd.concat([joined_df, df], axis=1, join="outer")
    return joined_df


def process_multiqc_dataframe(df):
    """
    Process the MultiQC dataframe by splitting the values in the first column by "." and setting it as the index.
    This function is required to handle the output files form quast & checkv to bring them to the same standard as MultiQC.

    Args:
        df (pd.DataFrame): The MultiQC dataframe to be processed.

    Returns:
        pd.DataFrame: The processed MultiQC dataframe.
    """
    # check if '.' are present in the first columns
    first_element = df[df.columns[0]].iloc[0]
    if "." in first_element and not re.search(r"\.\d", first_element):
        # split the first column by '.' and take the first part
        df[df.columns[0]] = df[df.columns[0]].str.split(".").str[0]
    df.set_index(df.columns[0], inplace=True)
    return df


def process_failed_contig_dataframe(df):
    """
    Process the failed contig dataframe by converting columns to string type,
    splitting values in "Cluster" and "Iteration" columns, creating a new "id" column,
    setting "id" as the index, and returning the index series.

    Args:
        df (pd.DataFrame): The dataframe containing the failed contig data.

    Returns:
        pandas.Index: The index series of the processed dataframe.
    """
    df = df.astype(str)
    df["id"] = df["sample name"] + "_" + df["cluster"] + "_" + df["step"]
    df.set_index("id", inplace=True)
    return df


def filter_files_of_interest(multiqc_data, files_of_interest):
    """
    Filters the multiqc_data list to include only files whose stem contains any of the files_of_interest.

    Args:
        multiqc_data (list): List of file paths.
        files_of_interest (list): List of file names to filter by.

    Returns:
        list: Filtered list of file paths.
    """
    file_list = [file for file in multiqc_data if files_of_interest in file.stem]
    if file_list:
        logger.debug("Files of interest found: %s", file_list)
    if len(file_list) > 1:
        logger.warning(
            "Multiple files of interest were found: %s for %s",
            file_list,
            files_of_interest,
        )
        logger.warning("Taking the first one: %s", file_list[0])
        return file_list[0]
    if len(file_list) == 0:
        return []
    return file_list[0]


def read_data(directory, file_columns, process_dataframe):
    """
    This function reads data from multiple files and processes it.

    Args:
    directory: The directory where the files are located.
    files_of_interest: A list of filenames that we are interested in.
    process_dataframe: A function that processes a dataframe.

    Returns:
    A dataframe that contains the processed data from all the files of interest.
    """
    logger.info("Reading data from %s", directory)
    multiqc_data = [file for file in directory.glob("multiqc_*.txt")]

    multiqc_samples_df = pd.DataFrame()
    for file_name, column_names in file_columns.items():
        files_of_interest = filter_files_of_interest(multiqc_data, file_name)
        if not files_of_interest:
            logger.info("No files of interest were found for %s in %s", file_name, directory)
            continue

        df = read_file_to_dataframe(files_of_interest)
        if df.empty:
            logger.warning("The file %s was empty!", files_of_interest)
            continue
        df = process_dataframe(df)
        df = filter_and_rename_columns(df, column_names)
        multiqc_samples_df = join_dataframe(multiqc_samples_df, [df])

    return multiqc_samples_df


def get_header(comment_dir, header_file_name):
    """
    Get the header from a header file located in the specified comment directory.

    Args:
        comment_dir (str): The directory where the header file is located.
        header_file_name (str): The name of the header file.

    Returns:
        list: The header as a list of strings.
    """
    header = []
    if comment_dir:
        header_file_path = f"{comment_dir}/{header_file_name}"
        if check_file_exists(header_file_path):
            header = read_header_file(header_file_path)
    return header


def dynamic_split(row, delimiter="_"):
    """
    Splits a string into multiple parts based on the specified delimiter.

    Args:
        row (str): The string to be split.
        delimiter (str, optional): The delimiter to use for splitting the string. Defaults to "_".

    Returns:
        list: A list of the split parts.
    """
    parts = row.split(delimiter)
    return parts


def compute_quast_metrics(quast_df):
    """
    Compute additional metrics based on QUAST output and add them to the DataFrame.

    Args:
        quast_df (pd.DataFrame): The QUAST output DataFrame.

    Returns:
        pd.DataFrame: The updated DataFrame with additional computed metrics.
    """
    if quast_df.empty:
        return quast_df

    try:
        # Compute the number of N's based on the "# N's per 100 kbp" and "Largest contig" columns
        quast_df["(quast) # N's"] = (
            pd.to_numeric(quast_df["(quast) # N's per 100 kbp"]) * pd.to_numeric(quast_df["(quast) Largest contig"]) / 100000
        )
        quast_df["(quast) # N's"] = quast_df["(quast) # N's"].astype(int)

        # Compute the percentage of N's based on the "# N's per 100 kbp" column
        quast_df["(quast) % N's"] = round(pd.to_numeric(quast_df["(quast) # N's per 100 kbp"]) / 1000, 2)

        # Keep only the relevant columns
        quast_df = quast_df[["(quast) # N's", "(quast) % N's", "(quast) # N's per 100 kbp"]]
    except KeyError as e:
        logger.warning(f"Missing column in QUAST output: {e}")
    except Exception as e:
        logger.error(f"Error computing QUAST metrics: {e}")

    return quast_df


def process_blast_dataframe(blast_df, blast_header=None, output_file=None):
    """
    Process the BLAST output DataFrame.

    Args:
        blast_df (pd.DataFrame): The BLAST output DataFrame.
        blast_header (list): A list of strings representing the header for the output file.
        output_file (str): The path to the output file.

    Returns:
        pd.DataFrame: The processed BLAST DataFrame.
    """
    if blast_df.empty:
        logger.warning("The BLAST DataFrame is empty.")
        return blast_df

    try:
        # Set the column names
        blast_df.columns = BLAST_COLUMNS

        # Filter for the best hit per contig and keep only the best hit
        blast_df = blast_df.sort_values("bitscore", ascending=False).drop_duplicates("query")

        # Process the DataFrame
        blast_df = generate_indexed_df(blast_df, "blast", "query", blast_header, output_file)

        # Make everything a string for the annotation
        blast_df = blast_df.astype(str)
    except Exception as e:
        logger.error(f"Error processing BLAST DataFrame: {e}")

    return blast_df


def process_annotation_dataframe(annotation_df, blast_header=None, output_file=None):
    """
    Process the annotation DataFrame.

    Args:
        annotation_df (pd.DataFrame): The annotation DataFrame.
        blast_header (list): A list of strings representing the header for the output file.
        output_file (str): The path to the output file.

    Returns:
        pd.DataFrame: The processed annotation DataFrame.
    """
    if annotation_df.empty:
        logger.warning("The annotation DataFrame is empty.")
        return annotation_df

    try:
        # Set the column names
        annotation_df.columns = BLAST_COLUMNS

        # Filter for the best hit per contig and keep only the best hit
        annotation_df = annotation_df.sort_values("bitscore", ascending=False).drop_duplicates("query")

        # Extract all key-value pairs into separate columns
        annotation_df = extract_annotation_data(annotation_df)

        # Remove subject title:
        annotation_df.drop(columns=["subject title"], inplace=True)

        annotation_df["% contig aligned"] = round((annotation_df["length"] / annotation_df["qlen"]) * 100, 2)

        # Process the DataFrame
        annotation_df = generate_indexed_df(annotation_df, "annotation", "query", blast_header, output_file)

        # Make everything a string for the annotation
        annotation_df = annotation_df.astype(str)
    except Exception as e:
        logger.error(f"Error processing annotation DataFrame: {e}")

    return annotation_df


def extract_annotation_data(df):
    # Extract all key-value pairs into separate columns
    df_extracted = df["subject title"].apply(parse_annotation_data).apply(pd.Series)
    return pd.concat([df, df_extracted], axis=1)


def parse_annotation_data(annotation_str):
    annotation_dict = {}
    pattern = r'(?P<key>\w+)\s*=\s*"?([^";]+)"?'
    matches = re.findall(pattern, annotation_str)
    for key, value in matches:
        annotation_dict[key] = value
    return annotation_dict


def generate_df(table_files, header_name=False, output=False, **kwargs):
    """
    Handle multiple table files and perform concatenation and writing to output file if specified.

    Args:
        table_files (list): List of table file paths.
        header_name (bool, optional): Flag to include header name in the output file. Defaults to False.
        output (str, optional): Output file path. Defaults to False.
        **kwargs: Additional keyword arguments for concatenation.

    Returns:
        pd.DataFrame: Concatenated table data.

    """
    result_df = pd.DataFrame()
    if table_files:
        result_df = concat_table_files(table_files, **kwargs)
    if output:
        write_dataframe(result_df, output, header_name)
    return result_df


def generate_indexed_df(df: pd.DataFrame, prefix: str, column_to_split: str, header=False, output=False) -> pd.DataFrame:
    """
    Handle the given dataframe by adding a prefix to column names, splitting a specific column,
    and generating an ID based on sample, cluster, and step information.

    Args:
        df (pd.DataFrame): The input dataframe.
        prefix (str): The prefix to add to column names.
        column_to_split (str): The column to split.
        header (bool, optional): Whether the dataframe has a header. Defaults to False.
        output (bool, optional): Whether to write the resulting dataframe to a file. Defaults to False.

    Returns:
        pd.DataFrame: The processed dataframe.
    """
    result_df = pd.DataFrame()
    if not df.empty:
        logger.info("Handling dataframe: %s", prefix)
        df = df.add_prefix(f"({prefix}) ")

        # Apply the dynamic split function to each row in the column
        split_data = df[f"({prefix}) {column_to_split}"].apply(dynamic_split).apply(pd.Series)
        # take the first three columns & rename
        split_data = split_data.iloc[:, :3]
        try:
            split_data.rename(
                columns={
                    split_data.columns[0]: "sample",
                    split_data.columns[1]: "cluster",
                    split_data.columns[2]: "step",
                },
                inplace=True,
            )
        except IndexError as e:
            logger.warning(
                "Unable to split up the file %s, at column %s \n ERROR: %s",
                prefix,
                column_to_split,
                e,
            )
            return result_df
        df = pd.concat([df, split_data], axis=1)
        df["step"] = df["step"].str.split(".").str[0]
        df["id"] = df["sample"] + "_" + df["cluster"] + "_" + df["step"]
        df = df.set_index("id")
        result_df = df
    if output:
        write_dataframe(result_df, output, header)
    result_df.drop(
        columns=["sample", "cluster", "step", f"({prefix}) {column_to_split}"],
        inplace=True,
    )
    return result_df


def filter_constrain(df, column, value):
    """
    Filter a dataframe based on a column and a regex value.

    Args:
        df (pd.DataFrame): The dataframe to be filtered.
        column (str): The column to filter on.
        regex_value (str): The regex value to filter on.

    Returns:
        pd.DataFrame, pd.DataFrame: The filtered dataframe with the regex value and the filtered dataframe without the regex value.
    """
    # Find rows with the regex value
    locations = df[column].str.contains(value) | df["step"].str.contains("constrain")

    # Filter
    df_with_value = df[locations]
    df_without_value = df[~locations]
    # Remove from column
    df_with_value.loc[:, column] = df_with_value[column].str.replace(value, "")
    df_with_value.loc[:, "index"] = df_with_value["index"].str.replace(value, "")
    return df_without_value, df_with_value


def drop_columns(df, columns):
    """
    Try to drop columns from a dataframe and return the dataframe.

    Args:
        df (pd.DataFrame): The dataframe to drop columns from.
        columns (list): The list of columns to drop.

    Returns:
        pd.DataFrame: The dataframe with the dropped columns.
    """
    result = df.drop(columns=[column for column in columns if column in df.columns])
    return result.copy()


def split_index_column(df):
    """
    Split the index column of the DataFrame into separate columns for sample name, cluster, and step.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The updated DataFrame with separate columns for sample name, cluster, and step.
    """
    df_copy = df.copy()
    df_copy = df_copy.reset_index().rename(columns={df_copy.index.name: "index"})
    df_copy = df_copy[df_copy["index"].str.contains("_", na=False)]
    df_copy[["sample name", "cluster", "step"]] = df_copy["index"].str.split("_", n=3, expand=True)
    return df_copy


def reorder_columns(df, columns):
    """
    Try to reorder columns in a dataframe and return the dataframe.

    Args:
        df (pd.DataFrame): The dataframe to reorder columns in.
        columns (list): The list of columns to reorder.

    Returns:
        pd.DataFrame: The dataframe with the reordered columns.
    """
    df = df[[column for column in columns if column in df.columns] + df.columns.difference(columns, sort=False).tolist()]
    return df


def reorder_rows(dataframe):
    """
    Reorder the rows in the DataFrame based on the ranking of the steps.

    Args:
        dataframe (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The reordered DataFrame.
    """

    df = dataframe.copy()
    ordered_list = ["constrain"] + [f"it{i}" for i in range(100, 0, -1)] + ["itvariant-calling", "consensus", "singleton"]
    rank_dict = {step: rank for rank, step in enumerate(ordered_list, start=1)}

    # Sort the DataFrame by 'step' based on the ranking dictionary
    df["rank"] = df["step"].map(rank_dict)
    df = df.sort_values(["sample name", "cluster", "rank"])

    return df


def filter_contigs(dataframe, level="normal") -> pd.DataFrame:
    """
    Filter contigs for each sample to only include those:
        - latest step of each cluster
        - Only annotated ones (if annotation is available)

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    df = reorder_rows(dataframe)

    # Select the first occurrence of each group
    df_snip = df.groupby(["sample name", "cluster"]).head(1).reset_index(drop=True)
    if not "annotation" in df.columns:
        logger.debug("Removed %d rows", len(df.index) - len(df_snip.index))
        return df_snip

    # Filter for annotated contigs
    df_snip = df_snip[df_snip["annotation"].notnull()]

    if level != "strict":
        logger.debug("Removed %d rows", len(df.index) - len(df_snip.index))
        return df_snip

    if [
        "(annotation) % contig aligned",
        "(annotation) bitscore",
    ] not in df_snip.columns:
        logger.warning("Columns for filtering are not present in the dataframe")
        return df_snip

    df_snip["sort"] = df_snip["(annotation) % contig aligned"] * df_snip["(annotation) bitscore"]
    df_snip = df_snip.sort_values(["sample name", "cluster", "sort"], ascending=[True, True, False])

    # Filter for the latest step of each cluster
    group_cols = [col for col in df_snip.columns if col in GROUPING_COLUMNS]
    logger.debug("Grouping columns for filtering: %s", group_cols)
    df_snip = df_snip.groupby([["sample name" + group_cols]]).tail(1).reset_index(drop=True)

    logger.debug("Removed %d rows", len(df.index) - len(df_snip.index))
    return df_snip


def coalesce_constrain(dataframe):
    """
    Fill missing values in the dataframe based on the group values.

    Args:
        dataframe (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The DataFrame with filled missing values.
    """

    df = reorder_rows(dataframe)
    grouping_cols = ["sample name", "cluster"]
    coalesce_columns = df.columns.difference(grouping_cols)
    result = df.copy()
    result[coalesce_columns] = result.groupby(grouping_cols)[coalesce_columns].transform(fill_group_na)

    return result.query('step == "constrain"')


def fill_group_na(s):
    return s.infer_objects(copy=False).ffill().bfill()


def read_bed(bed_files, selection):
    """
    Read bed files and concatenate them into a single dataframe.

    Args:
        bed_files (list): List of bed file paths.
        selection (list): List of approved bed_files.

    Returns:
        dict: {index: pd.Dataframe} dictionary of bed dataframes.
    """
    bed_data = {}
    for sample in selection:
        bed_sel = [bed_file for bed_file in bed_files if sample_in_file(sample, bed_file)]
        logger.debug("Bed files for %s: %s", sample, bed_sel)
        if not bed_sel:
            logger.warning("No bed files were found for %s", sample)
            continue
        bed_file = bed_sel[0]
        if check_file_exists(bed_file):
            logger.debug("Converting bed file %s to coverage", bed_file)
            df = pd.read_csv(
                bed_file,
                sep="\t",
                compression="gzip",
                names=["chromosome", "position", "end", "coverage"],
            )
            bed_data[sample] = df[["position", "coverage"]].copy()

    if len(bed_data) != len(selection):
        logger.warning("Not all bed files were read!")
        logger.warning("Missing bed files: %s", set(selection) - set(bed_data.keys()))
    return bed_data


def sample_in_file(sample, file_name) -> bool:
    """
    Check if the sample name is in the file name.

    Args:
        sample (str): The sample name.
        file (str): The file name.

    Returns:
        bool: True if the sample name is in the file name, False otherwise.
    """
    sample_upper = sample.upper()
    file_upper = str(file_name).upper()
    logger.debug("sample: %s - file: %s", sample_upper, file_upper)

    return sample_upper in file_upper or sample_upper in file_upper.replace("_", "-") or sample_upper in file_upper.replace("-", "_")


def remove_keys(d, keys):
    """
    Remove keys from a dictionary.

    Args:
        d (dict): The dictionary to remove keys from.
        keys (list): The list of keys to remove.

    Returns:
        dict: The dictionary with the specified keys removed.
    """
    return {k: v for k, v in d.items() if k not in keys}


def create_constrain_summary(df_constrain: pd.DataFrame, file_columns: List[Union[str, Dict[str, str]]]) -> pd.DataFrame:
    """
    Create a summary table for the constrain data.

    Args:
        df_constrain (pd.DataFrame): The constrain data DataFrame.
        file_columns (List): A columns dictionary with old names & new names.

    Returns:
        pd.DataFrame: The constrain summary table.
    """

    # Filter only for columns of interest
    # Some columns were already renamed, so we get the new values of them based on the original naming of mqc
    dic_columns = {}
    for item in file_columns:
        if isinstance(item, dict):
            dic_columns.update(item)
        else:
            dic_columns[item] = item

    columns_of_interest = [dic_columns[key] for key in CONSTRAIN_GENERAL_STATS_COLUMNS if key in dic_columns.keys()]

    if not columns_of_interest:
        logger.warning("No columns of interest were found to create the constrain summary table!")
        return pd.DataFrame()

    columns_of_interest = [
        "sample name",
        "species",
        "segment",
        "cluster",
        "definition",
        "qlen",  # length of the query sequence will have to be renamed
    ] + columns_of_interest

    df_columns = df_constrain.columns.tolist()

    present_columns = []
    for name in columns_of_interest:
        if name in df_columns:  # Check for an exact match first
            present_columns.append(name)
        else:  # If no exact match, try approximate match
            matches = [col for col in df_columns if name in col]
            if matches:
                matched_column = matches[0]
                present_columns.append(matched_column)

    df_constrain = df_constrain[present_columns]

    if df_constrain.empty:
        return df_constrain

    if "(blast) qlen" in df_constrain.columns:
        df_constrain = df_constrain.rename(columns={"(blast) qlen": "consensus length"})

    # Reformat dataframe to long based on following:
    #   Species & Segment
    #   Species
    #   ID (Cluster)
    df_constrain.loc[:, "idgroup"] = df_constrain.apply(
        lambda row: (
            f"{row['species']} ({row['segment']})"
            if "segment" in df_constrain.columns and pd.notnull(row["species"]) and pd.notnull(row["segment"])
            else (row["species"] if "species" in df_constrain.columns and pd.notnull(row["species"]) else row["cluster"])
        ),
        axis=1,
    )
    df_constrain = df_constrain.rename(columns={"cluster": "Constrain id"})

    # Remove columns that are not needed anymore
    df_constrain = drop_columns(df_constrain, ["species", "segment"])

    # Convert dataframe to long and then extra wide
    df_long = df_constrain.melt(id_vars=["idgroup", "sample name"], var_name="variable", value_name="Value")
    # Remove rows with NaN values & duplicates
    df_long = df_long.dropna()
    df_long = df_long.drop_duplicates()
    df_long["grouped variable"] = df_long["idgroup"] + " - " + df_long["variable"]
    df_long.drop(columns=["idgroup", "variable"], inplace=True)
    # Convert to wide format
    df_wide = df_long.pivot(index=["sample name"], columns="grouped variable", values="Value")
    df_wide.reset_index(inplace=True)

    return df_wide


def load_custom_data(args) -> List[pd.DataFrame]:
    """
    Load custom data from files and process it to a list of dataframes.
    """

    # Clusters overview - mini multiqc module
    if args.clusters_summary:
        clusters_df = generate_df(args.clusters_summary)
        if not clusters_df.empty:
            clusters_df.set_index("Sample name", inplace=True)
            module = mqc.BaseMultiqcModule(name="Cluster Summary", anchor=Anchor("cluster-summary"))
            module.general_stats_addcols(clusters_df.to_dict(orient="index"))

            # Custom barplot -  Clusters sample
            plot_df = reorder_columns(clusters_df.copy(), ["# Clusters", "Filtered # clusters", "Total # clusteres"])
            plot = bargraph.plot(data=plot_df.to_dict(orient="index"), pconfig=CLUSTER_PCONFIG)
            module.add_section(
                anchor=Anchor("cluster-summary"), plot=plot, description="Number of identified contig clusters per sample after assembly."
            )
            mqc.report.modules.append(module)

    # General Stats - Sample metadata
    if args.sample_metadata:
        metadata_df = generate_df([args.sample_metadata])
        if not metadata_df.empty:
            sample_col = [col for col in metadata_df.columns if "sample" in col.lower()][0]
            metadata_df.set_index(sample_col, inplace=True)
            module = mqc.BaseMultiqcModule(name="Sample metadata", anchor=Anchor("custom_data"))
            content = metadata_df.to_dict(orient="index")
            module.general_stats_addcols(content)
            mqc.report.modules.append(module)

    # CLuster table - Checkv summary
    checkv_df = generate_df(args.checkv_files)
    if not checkv_df.empty:
        checkv_df = generate_indexed_df(checkv_df, "checkv", "contig_id")

    # Cluster table - Quast summary
    quast_df = read_in_quast(args.quast_files)
    if not quast_df.empty:
        quast_df = generate_indexed_df(quast_df, "quast", "Assembly")

        # Most of the columns are not good for a single contig evaluation
        quast_df["(quast) # N's"] = (
            pd.to_numeric(quast_df["(quast) # N's per 100 kbp"]) * pd.to_numeric(quast_df["(quast) Largest contig"]) / 100000
        )
        quast_df = quast_df.astype({"(quast) # N's": int})
        quast_df["(quast) % N's"] = round(pd.to_numeric(quast_df["(quast) # N's per 100 kbp"]) / 1000, 2)
        quast_df = quast_df[["(quast) # N's", "(quast) % N's", "(quast) # N's per 100 kbp"]]

    # Cluster table - Blast summary
    blast_df = generate_df(args.blast_files, header=None)
    if not blast_df.empty:
        blast_df = process_blast_dataframe(blast_df)

    # Cluster table - mmseqs easysearch summary (annotation section)
    annotation_df = generate_df(args.annotation_files, header=None)
    if not annotation_df.empty:
        annotation_df = process_annotation_dataframe(annotation_df)

    return [checkv_df, quast_df, blast_df, annotation_df]


def get_module_data(mqc: object, module: str) -> Dict[str, any]:
    """
    Attempt to get data for a module that might be a partial match.

    Args:
        mqc (object): MultiQC object.
        module (str): Module name to search for.

    Returns:
        Dict[str, Any]: Module data if found, otherwise an empty dict.
    """
    if module in mqc.list_modules():
        data = mqc.get_module_data(module)
        return data

    module_basename = module.split("_")[0]
    if module_basename in mqc.list_modules():
        all_data = mqc.get_module_data(module_basename)
        if all_data:
            matching_key = next((key for key in all_data.keys() if module in key), None)
            if matching_key:
                logger.debug("Data found for %s in MultiQC under key %s", module, matching_key)
                return all_data[matching_key]
    return {}


def flatten(xss):
    return [x for xs in xss for x in xs]


def handle_module_data(
    mqc: object,
    module: str,
    section: Union[
        str,
        List[
            Union[
                str,
                Dict[str, str],
                Dict[
                    str,
                    List[Union[str, Dict[str, str]]],  # module has multiple sections
                ],
            ]
        ],
    ],
) -> Tuple[list[pd.DataFrame], List[Union[str, Dict[str, str]]]]:
    def check_section_exists(module_data: Dict, section_key: str) -> bool:
        """Check if a section exists in the module data."""
        return any(section_key in key for key in module_data.keys())

    # Get all module data first
    all_module_data = mqc.get_module_data(module=module)
    logger.debug("All data for %s, %s", module, all_module_data)

    if not all_module_data:
        logger.warning("No data found for module: %s", module)
        return [pd.DataFrame()], []

    # Empty section, so we take all values from module
    if isinstance(section, str):
        if not section:
            return [pd.DataFrame.from_dict(all_module_data, orient="index")], []
        else:
            if check_section_exists(all_module_data, section):
                return [pd.DataFrame.from_dict(all_module_data[section], orient="index")], []
            else:
                logger.warning("Section %s not found in module %s", section, module)
                return [pd.DataFrame()], []
    elif isinstance(section, list):
        if isinstance(section[0], str):
            # Section refers to column names
            return [filter_and_rename_columns(pd.DataFrame.from_dict(all_module_data, orient="index"), section)], section
        elif isinstance(section[0], dict):
            first_value = next(iter(section[0].values()))
            if isinstance(first_value, str):
                # Already at column level
                return [filter_and_rename_columns(pd.DataFrame.from_dict(all_module_data, orient="index"), section)], section
            elif isinstance(first_value, list):
                # Section could have multiple sections:
                result_df = []
                result_list = []
                for subsection in section:
                    subsection_data, columns = handle_module_data(mqc, module, subsection)
                    if isinstance(subsection_data, pd.DataFrame):
                        result_df.extend(subsection_data)
                    if isinstance(columns, list):
                        result_list.extend(columns)
                return result_df, result_list
    elif isinstance(section, dict):
        # We just have {section: List[Union[colname, Dict[colname, colrename]]}
        section_name, columns = next(iter(section.items()))
        if check_section_exists(all_module_data, section_name):
            data = pd.DataFrame.from_dict(all_module_data[section_name], orient="index")
            return [filter_and_rename_columns(data, columns)], columns
        else:
            logger.warning("Section '%s' not found in module '%s'", section_name, module)
            return [pd.DataFrame()], []

    logger.warning(
        "Unsupported section type from module %s: %s for %s",
        module,
        type(section),
        section,
    )
    return [pd.DataFrame()], []


def handle_general_stats(columns: List[Union[str, Dict[str, str]]]) -> tuple[pd.DataFrame, List]:
    """
    Handle general stats data from MultiQC.

    Args:
        columns (List[Union[str, Dict[str,str]]): List of columns to filter and rename.

    Returns:
        pd.DataFrame: The filtered and renamed DataFrame.
    """

    df = pd.DataFrame.from_dict(mqc.get_general_stats_data(), orient="index")

    return [filter_and_rename_columns(df, columns)], []


def extract_mqc_data(table_headers: Union[str, Path]) -> Optional[pd.DataFrame]:
    """
    Extract data from MultiQC output files.

    Args:
        mqc (object): MultiQC object.
        table_headers (Union[str, Path]): Path to the table headers file.

    Returns:
        pd.DataFrame: Extracted data
    """
    result = pd.DataFrame()
    module_selection = get_module_selection(table_headers)
    av_modules = mqc.list_modules()
    data = []
    columns_result = []

    for module, section in module_selection.items():
        if module == "general_stats":
            logger.info("Extracting general stats data from multiqc")
            module_data, columns = handle_general_stats(section)

        elif module not in av_modules:
            logger.warning("Module %s is not available in MultiQC, skipping extraction", module)
            continue

        else:
            logger.info("Extracting %s data from multiqc", module)
            module_data, columns = handle_module_data(mqc, module, section)
            logger.debug("Data for %s: %s", module, module_data)

        data.extend(module_data)
        columns_result.extend(columns)

    logger.debug("Data list: %s", data)

    return join_dataframe(result, data) if data else result, columns_result


def reformat_custom_df(df):
    """
    Reformat the custom dataframe.
    """
    # Keep only those rows we can split up in sample, cluster, step
    logger.info("Splitting up the index column in sample name, cluster, step")
    df = split_index_column(df)

    # Reorder the columns
    logger.info("Reordering columns")
    final_columns = ["index", "sample name", "cluster", "step"] + [
        column
        for group in [
            "annotation",
            "mas-screen",
            "blast",
            "checkv",
            "QC check",
            "quast",
        ]
        for column in df.columns
        if group in column
    ]
    return reorder_columns(df, final_columns)


def reformat_constrain_df(df, file_columns, args):
    """
    Reformat the constrain dataframe.
    """
    # Separate table for mapping constrains
    if df.empty:
        return mqc, df

    # Add constrain metadata to the mapping constrain table
    constrain_meta = generate_df([args.mapping_constrains])

    # drop unwanted columns & reorder
    constrain_meta = drop_columns(constrain_meta, ["sequence", "samples"])
    df = df.merge(constrain_meta, how="left", left_on="cluster", right_on="id")
    df = reorder_columns(
        df,
        [
            "index",
            "sample name",
            "cluster",
            "step",
            "species",
            "segment",
            "definition",
        ],
    )

    # add mapping summary to sample overview table in ... wide format with species & segment combination
    logger.info("Creating mapping constrain summary (wide) table")
    mapping_constrains_summary = create_constrain_summary(df, file_columns).set_index("sample name")
    write_dataframe(mapping_constrains_summary, "mapping_constrains_summary.tsv", [])
    if not mapping_constrains_summary.empty:
        # Add to mqc
        module = mqc.BaseMultiqcModule(name="Mapping Constrains Summary", anchor=Anchor("custom_data"))
        content = mapping_constrains_summary.to_dict(orient="index")
        module.general_stats_addcols(content)
        mqc.report.modules.append(module)

    logger.info("Coalescing columns")
    coalesced_constrains = coalesce_constrain(df)
    return coalesced_constrains


def write_results(contigs_mqc, constrains_mqc, args) -> int:
    """
    Write the results to files.
    """

    if not contigs_mqc.empty:
        logger.info("Writing Unfiltered Denovo constructs table file: contigs_all.tsv")
        write_dataframe(contigs_mqc, "contigs_all.tsv", [])
        contigs_mqc.set_index("index", inplace=True)
        mqc.add_custom_content_section(
            name="Denovo Construct Overview",
            anchor=Anchor("contigs_all"),
            description="The table below shows the overview of the denovo constructs with refinement.",
            plot=table.plot(data=contigs_mqc.to_dict(orient="index")),
        )

    if not constrains_mqc.empty:
        logger.info("Writing Unfiltered Mapping constructs table file: mapping_all.tsv")
        write_dataframe(constrains_mqc, "mapping_all.tsv", [])
        constrains_mqc.set_index("index", inplace=True)
        mqc.add_custom_content_section(
            name="Mapping Construct Overview",
            anchor=Anchor("mapping_all"),
            description="The table below shows the overview of the mapping constructs with refinement.",
            plot=table.plot(data=constrains_mqc.to_dict(orient="index")),
        )

    # TODO correctly insert metadata of:
    #   -  versions
    #   -  citations
    #   -  parameters
    #   -  methods_description
    # TODO: Include only final iteration in plots
    mqc.write_report(make_data_dir=True, data_format="tsv", export_plots=False)

    return 0


def generate_ignore_sample_pattern(df: pd.DataFrame) -> pd.DataFrame:
    """ """
    base_pattern = r".*_"
    ignores = ["consensus", "singleton"]

    it_patterns = [0]
    # Find all 'it' patterns in index
    for idx in df.index:
        parts = idx.split("_")
        if parts and parts[-1].startswith("it") and parts[-1][2:].isdigit():
            it_patterns.append(int(parts[-1][2:]))

    max_number = max(it_patterns)

    if max_number == 0:
        return f'{base_pattern}({"|".join(ignores)})$'

    ignores.append("itvariant_calling")
    ignores.append([f"it{i}" for i in range(0, max_number)])

    return f'{base_pattern}({"|".join(ignores)})$'


def main(argv=None):
    """
    Main function for creating custom tables for MultiQC.

    Args:
        argv (list): List of command line arguments.

    Returns:
        int: Exit code.
    """
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level,
        format="[%(asctime)s - %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # 1. Run MQC with correct config
    mqc.parse_logs(
        args.multiqc_files,
        args.multiqc_config,
    )
    for module in mqc.list_modules():
        module_data = mqc.get_module_data(module)
        logger.info("Data for %s: %s", module, module_data.keys())

    # 2. Extract MQC data
    mqc_custom_df, renamed_columns = extract_mqc_data(args.table_headers)
    if mqc_custom_df.empty:
        logger.warning("No data was found from MULTIQC to create the contig overview table! - Exiting")
        return 0

    # 3. Reset multiqc and rerun while removing iteration data.
    mqc.reset()
    mqc.parse_logs(args.multiqc_files, args.multiqc_config, ignore_samples=generate_ignore_sample_pattern(mqc_custom_df))

    # 2. Parse our custom files into the correct tables
    custom_tables = load_custom_data(args)

    # 3. Make our own summary excel
    # 3.1 Extract the MQC data

    # 3.2 Join with the custom contig tables
    mqc_custom_df = join_dataframe(mqc_custom_df, custom_tables)

    if mqc_custom_df.empty:
        logger.warning("No data was found to create the contig overview table!")
        return 0

    # 3.3 reformat the dataframe
    mqc_custom_df = reformat_custom_df(mqc_custom_df)

    # 3.4 split up denovo constructs and mapping (-CONSTRAIN) results
    logger.info("Splitting up denovo constructs and mapping (-CONSTRAIN) results")
    contigs_mqc, constrains_mqc = filter_constrain(mqc_custom_df, "cluster", "-CONSTRAIN")

    coalesced_constrains = reformat_constrain_df(constrains_mqc, renamed_columns, args)

    write_results(contigs_mqc, coalesced_constrains, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
