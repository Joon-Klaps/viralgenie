#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import argparse
import csv
import logging
import json
import os
import re
import sys
from constant_variables import (
    BLAST_COLUMNS,
    CONSTRAIN_GENERAL_STATS_COLUMNS,
    FILES_OF_INTEREST,
    CLUSTER_PCONFIG,
)
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


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Provide a command line tool to combine individual log & summary files which we will pass down to multiqc subsequently.",
        epilog="Example: python custom_multiqc.py --clusters_summary file1,file2,file3,... ",
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


def dynamic_split(index_str):
    """
    Split the index string into sample name, cluster, and step.
    """
    parts = index_str.split("_")
    return parts[0], parts[1], "_".join(parts[2:])


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
        blast_df = generate_indexed_df(blast_df, "blast", "query")

        # Make everything a string for the annotation
        blast_df = blast_df.astype(str)
    except Exception as e:
        logger.error(f"Error processing BLAST DataFrame: {e}")

    return blast_df


def process_annotation_dataframe(annotation_df):
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
        annotation_df = generate_indexed_df(annotation_df, "annotation", "query")

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


def generate_indexed_df(df: pd.DataFrame, prefix: str = None, column_to_split: str = "index") -> pd.DataFrame:
    """
    Handle the given dataframe by adding a prefix to column names, splitting a specific column,
    and generating an ID based on sample, cluster, and step information.

    Args:
        df (pd.DataFrame): The input dataframe.
        prefix (str): The prefix to add to column names.
        column_to_split (str): The column to split.

    Returns:
        pd.DataFrame: The processed dataframe.
    """
    result_df = pd.DataFrame()
    if df.empty:
        return result_df

    split_column = column_to_split
    if prefix:
        logger.info("Handling dataframe: %s", prefix)
        df = df.add_prefix(f"({prefix}) ")
        split_column = f"({prefix}) {column_to_split}"

    df = split_index_column(df, prefix, split_column)

    df["step"] = df["step"].str.split(".").str[0]
    df["id"] = df["sample"] + "_" + df["cluster"] + "_" + df["step"]
    df = df.set_index("id")
    result_df = df.copy()
    result_df.drop(
        columns=["sample", "cluster", "step", split_column],
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


def drop_columns(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
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


def split_index_column(df: pd.DataFrame, prefix: str = None, split_column: str = "index") -> pd.DataFrame:
    """
    Split the index column of the DataFrame into separate columns for sample name, cluster, and step.

    Args:
        df (pd.DataFrame): The input DataFrame.
        prefix (str, optional): A prefix for the DataFrame.
        split_column (str, optional): The column to split. Defaults to "index".

    Returns:
        pd.DataFrame: The updated DataFrame with separate columns for sample name, cluster, and step.
    """
    df_copy = df.copy()
    # Reset the index and rename the index column
    df_copy = df_copy.reset_index(drop=True).rename(columns={df_copy.index.name: split_column})
    df_copy = df_copy[df_copy[split_column].str.contains("_", na=False)]

    # Apply the dynamic split function to each row in the column
    split_data = df_copy[split_column].apply(dynamic_split).apply(pd.Series)

    # Take the first three columns and rename them
    split_data = split_data.iloc[:, :3]
    split_data.rename(
        columns={
            split_data.columns[0]: "sample",
            split_data.columns[1]: "cluster",
            split_data.columns[2]: "step",
        },
        inplace=True,
    )

    df_copy = drop_columns(df_copy, ["sample", "cluster", "step"])
    # Concatenate the original DataFrame and the split data
    df_copy = pd.concat([df_copy, split_data], axis=1)

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
    df["rank"] = df["step"].replace(rank_dict)
    df = df.sort_values(["sample", "cluster", "rank"])

    return df


def coalesce_constrain(dataframe):
    """
    Fill missing values in the dataframe based on the group values.

    Args:
        dataframe (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The DataFrame with filled missing values.
    """

    df = reorder_rows(dataframe)
    grouping_cols = ["sample", "cluster"]
    coalesce_columns = df.columns.difference(grouping_cols)
    result = df.copy()
    result[coalesce_columns] = result.groupby(grouping_cols)[coalesce_columns].transform(fill_group_na)

    return result.query('step == "constrain"')


def fill_group_na(s):
    return s.infer_objects(copy=False).ffill().bfill()


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
        "sample",
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
    df_long = df_constrain.melt(id_vars=["idgroup", "sample"], var_name="variable", value_name="Value")
    # Remove rows with NaN values & duplicates
    df_long = df_long.dropna()
    df_long = df_long.drop_duplicates()
    df_long["grouped variable"] = df_long["idgroup"] + " - " + df_long["variable"]
    df_long.drop(columns=["idgroup", "variable"], inplace=True)
    # Convert to wide format
    df_wide = df_long.pivot(index=["sample"], columns="grouped variable", values="Value")
    df_wide.reset_index(inplace=True)

    return df_wide


def select_columns(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Try to select columns from a dataframe and return the dataframe.

    Args:
        df (pd.DataFrame): The dataframe to select columns from.
        columns (list): The list of columns to select.

    Returns:
        pd.DataFrame: The dataframe with the selected columns.
    """
    result = df[[column for column in columns if column in df.columns]]
    return result.copy()


def load_custom_data(args) -> List[pd.DataFrame]:
    """
    Load custom data from files and process it to a list of dataframes.
    """

    # Clusters overview - mini multiqc module
    if args.clusters_summary:
        clusters_df = filelist_to_df(args.clusters_summary)
        if not clusters_df.empty:
            clusters_df.set_index("Sample name", inplace=True)

            # Adding to general stats
            module = mqc.BaseMultiqcModule(name="Cluster Summary", anchor=Anchor("cluster-summary"))
            module.general_stats_addcols(clusters_df.to_dict(orient="index"))

            # Custom barplot -  Clusters sample
            plot_df = select_columns(clusters_df, ["# Clusters", "# Removed clusters"])
            plot = bargraph.plot(data=plot_df.to_dict(orient="index"), pconfig=CLUSTER_PCONFIG)
            module.add_section(
                anchor=Anchor("cluster-summary"), plot=plot, description="Number of identified contig clusters per sample after assembly."
            )
            mqc.report.modules.append(module)

    # General Stats - Sample metadata
    if args.sample_metadata:
        metadata_df = filelist_to_df([args.sample_metadata])
        if not metadata_df.empty:
            sample_col = [col for col in metadata_df.columns if "sample" in col.lower()][0]
            metadata_df.set_index(sample_col, inplace=True)
            module = mqc.BaseMultiqcModule(name="Sample metadata", anchor=Anchor("custom_data"))
            content = metadata_df.to_dict(orient="index")
            module.general_stats_addcols(content)
            mqc.report.modules.append(module)

    # CLuster table - Checkv summary
    checkv_df = filelist_to_df(args.checkv_files)
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
    blast_df = filelist_to_df(args.blast_files, header=None)
    if not blast_df.empty:
        blast_df = process_blast_dataframe(blast_df)

    # Cluster table - mmseqs easysearch summary (annotation section)
    annotation_df = filelist_to_df(args.annotation_files, header=None)
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
    """flatten a list"""
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
    """
    Extract data from multiqc modules based on a nested yml file structure for filtering.

    Args:
        mqc (object): MultiQC object.
        module (str): The module name.
        section (Union[str, List[Union[str, Dict[str, str]]]): The section to extract data from.

    Returns:
        Tuple[list[pd.DataFrame], List[Union[str, Dict[str, str]]]]: A list of dataframes and a list of columns.
    """

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
    df = drop_columns(df, ["index"])
    df["index"] = df.index

    df = split_index_column(df)

    # Reorder the columns
    logger.info("Reordering columns")
    final_columns = ["index", "sample", "cluster", "step"] + [
        column
        for group in [
            "annotation",
            "mash-screen",
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
    constrain_meta = filelist_to_df([args.mapping_constrains])

    # drop unwanted columns & reorder
    constrain_meta = drop_columns(constrain_meta, ["sequence", "samples"])
    df = df.merge(constrain_meta, how="left", left_on="cluster", right_on="id")
    df = reorder_columns(
        df,
        [
            "index",
            "sample",
            "cluster",
            "step",
            "species",
            "segment",
            "definition",
        ],
    )

    # add mapping summary to sample overview table in ... wide format with species & segment combination
    logger.info("Creating mapping constrain summary (wide) table")
    mapping_constrains_summary = create_constrain_summary(df, file_columns).set_index("sample")
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
        table_plot = contigs_mqc[~generate_ignore_samples(contigs_mqc)]
        mqc.add_custom_content_section(
            name="Denovo Construct Overview",
            anchor=Anchor("contigs_all"),
            description="The table below shows the overview of the denovo constructs with refinement.",
            plot=table.plot(data=table_plot.to_dict(orient="index")),
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
    #   -  Not all mapping data is in the general stats table, while it should be
    #   -  Double check for any other loss of information.
    mqc.write_report(make_data_dir=True, data_format="tsv", export_plots=False)

    return 0


def generate_ignore_samples(dataframe: pd.DataFrame) -> pd.Series:
    """
    Generate a Series of indices that are not part of the df_snip dataframe.

    Parameters:
    dataframe (pd.DataFrame): The input DataFrame to process.

    Returns:
    pd.Series: A Series containing the indices that are not in df_snip.
    """
    df = dataframe.copy()
    df = drop_columns(df, ["index"])
    df["index"] = df.index
    df = split_index_column(df)

    df = reorder_rows(df)

    # Filter for only the last iteration
    df_filter = df.groupby(["sample", "cluster"]).head(1).reset_index(drop=True)

    return df["index"][~df["index"].isin(df_filter["index"])]


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

    for module in [m for m in mqc.list_modules() if "viralgenie" not in m]:
        module_data = mqc.get_module_data(module)
        logger.info("Data for %s: %s", module, module_data.keys())

    # 2. Extract MQC data
    mqc_custom_df, renamed_columns = extract_mqc_data(args.table_headers)
    if mqc_custom_df.empty:
        logger.warning("No data was found from MULTIQC to create the contig overview table! - Exiting")
        return 0

    # 3. Reset multiqc and rerun while removing iteration data.
    mqc.reset()
    mqc.parse_logs(args.multiqc_files, args.multiqc_config, ignore_samples=generate_ignore_samples(mqc_custom_df))

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
