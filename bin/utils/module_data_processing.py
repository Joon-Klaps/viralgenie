#!/usr/bin/env python

import logging
import re
from typing import Dict, List, Union, Tuple, Optional, Any

import pandas as pd

from utils.constant_variables import BLAST_COLUMNS, CONSTRAIN_GENERAL_STATS_COLUMNS, COLUMN_MAPPING
from utils.file_tools import filelist_to_df
from utils.pandas_tools import (
    coalesce_constrain,
    drop_columns,
    generate_indexed_df,
    reorder_columns,
    reorder_rows,
    split_index_column,
    filter_and_rename_columns,
)

logger = logging.getLogger()
pd.options.mode.copy_on_write = True


def process_blast_df(blast_df):
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


def process_annotation_df(annotation_df):
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


def reformat_custom_df(df: pd.DataFrame, cluster_df: pd.DataFrame) -> pd.DataFrame:
    """
    Reformat the custom dataframe.
    """
    # Keep only those rows we can split up in sample, cluster, step
    logger.info("Splitting up the index column in sample name, cluster, step")
    df = drop_columns(df, ["index"])
    df["index"] = df.index

    df = split_index_column(df)

    if not cluster_df.empty:
        df = pd.merge(df, cluster_df, on=["sample", "cluster"], how="left")
        df.index = df["index"]

    # Reorder the columns
    logger.info("Reordering columns")
    final_columns = ["index", "sample", "cluster", "step"] + [
        column
        for group in [
            "annotation",
            "mash-screen",
            "blast",
            "checkv",
            "cluster",
            "quast",
        ]
        for column in df.columns
        if group in column
    ]
    return reorder_columns(df.dropna(subset=["step"]), list(dict.fromkeys(final_columns)))


def filter_constrain(dataframe, column, value):
    """
    Filter a dataframe based on a column and a regex value.

    Args:
        dataframe (pd.DataFrame): The dataframe to be filtered.
        column (str): The column to filter on.
        regex_value (str): The regex value to filter on.

    Returns:
        pd.DataFrame, pd.DataFrame: The filtered dataframe with the regex value and the filtered dataframe without the regex value.
    """
    df = dataframe.copy()
    # Find rows with the regex value
    locations = df[column].str.contains(value) | df["step"].str.contains("constrain")

    # Filter
    df_with_value = df[locations]
    df_without_value = df[~locations]
    # Remove from column
    df_with_value.loc[:, column] = df_with_value[column].str.replace(value, "")
    df_with_value.loc[:, "index"] = df_with_value["index"].str.replace(value, "")

    return df_without_value.dropna(axis=1, how="all"), df_with_value.dropna(axis=1, how="all")


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

    logger.debug("dic_columns: %s", dic_columns)

    columns_of_interest = [dic_columns.get(key, key) for key in CONSTRAIN_GENERAL_STATS_COLUMNS]

    logger.debug("columns_of_interest: %s", columns_of_interest)

    if not columns_of_interest:
        logger.warning("No columns of interest were found to create the constrain summary table!")
        return pd.DataFrame()

    columns_of_interest.extend(
        [
            "sample",
            "species",
            "segment",
            "cluster",
            "definition",
        ]
    )

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
        logger.warning("The constrain DataFrame is empty.")
        return df_constrain

    df_constrain = df_constrain.rename(columns=COLUMN_MAPPING)

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


def reformat_constrain_df(df, file_columns, args):
    """
    Reformat the constrain dataframe.
    """
    # Separate table for mapping constrains
    if df.empty:
        logger.warning("The constrain DataFrame is empty.")
        return df, df

    # Add constrain metadata to the mapping constrain table
    constrain_meta = filelist_to_df([args.mapping_constraints])

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
    mapping_constraints_summary = create_constrain_summary(df, file_columns).set_index("sample")

    logger.info("Coalescing columns")
    coalesced_constrains = coalesce_constrain(df)
    coalesced_constrains = drop_columns(coalesced_constrains, ["id", "selection", "rank"])
    return coalesced_constrains, mapping_constraints_summary


def generate_ignore_samples(dataframe: pd.DataFrame) -> pd.Series:
    """
    Generate a Series of indices that are not part of the df_snip dataframe.

    Parameters:
    dataframe (pd.DataFrame): The input DataFrame to process.

    Returns:
    pd.Series: A Series containing the indices that are not in df_snip.
    """
    df = dataframe.copy()
    df = split_index_column(df)

    df = reorder_rows(df)

    # Filter for only the last iteration
    df_filter = df.groupby(["sample", "cluster"]).head(1).reset_index(drop=True)

    return df["index"][~df["index"].isin(df_filter["index"])]


def add_prefix_to_values_dict(data: List[Union[str, Dict[str, str]]], prefix: str) -> List[Dict[str, str]]:
    updated_items = []
    for item in data:
        if isinstance(item, str):
            updated_items.append({item: f"({prefix}) {item}"})
        else:
            updated_items.extend({key: f"({prefix}) {value}"} for key, value in item.items())
    return updated_items


def check_section_exists(module_data: Dict, section_key: str) -> bool:
    """Check if a section exists in the module data."""
    return any(section_key in key for key in module_data.keys())


def extract_mqc_from_simple_section(all_module_data: Dict, section: Optional[str], module: str) -> Tuple[List[pd.DataFrame], List[Any]]:
    """Handle simple string or None section cases."""
    logger.debug("Extracting data from simple str %s", module)
    if not section:
        # Return all data if no specific section is specified
        return [pd.DataFrame.from_dict(all_module_data, orient="index")], []

    # Check if specific section exists
    if check_section_exists(all_module_data, section):
        return [pd.DataFrame.from_dict(all_module_data[section], orient="index")], []

    logger.warning(f"Section {section} not found in module {module}")
    return [pd.DataFrame()], []


def extract_mqc_from_list_section(all_module_data: Dict, section: List, module: str) -> Tuple[List[pd.DataFrame], List[Any]]:
    """Handle list-based section specifications."""
    logger.debug("Extracting data from list %s : %s", module, section)
    # Case for list of column names
    if all(not isinstance(item, dict) or not isinstance(list(item.values())[0], list) for item in section):
        full_df = pd.DataFrame.from_dict(all_module_data, orient="index")
        return [filter_and_rename_columns(full_df, section)], section

    # Handle nested section lists
    result_dfs = []
    result_columns = []
    for subsection in section:
        # Handle different types of subsections
        if isinstance(subsection, str):
            # Simple section name
            subsection_dfs, subsection_columns = extract_mqc_from_simple_section(all_module_data, subsection, module)
        if isinstance(subsection, list):
            # Simple section name
            subsection_dfs, subsection_columns = extract_mqc_from_list_section(all_module_data, subsection, module)
        elif isinstance(subsection, dict):
            # Dictionary-based section specification
            subsection_dfs, subsection_columns = extract_mqc_from_dict_section(all_module_data, subsection, module)
        else:
            # Unsupported subsection type
            logger.warning(f"Unsupported subsection type: {type(subsection)}")
            continue

        result_dfs.extend(subsection_dfs)
        result_columns.extend(subsection_columns)

    return result_dfs, result_columns


def extract_mqc_from_dict_section(all_module_data: Dict, section: Dict, module: str) -> Tuple[List[pd.DataFrame], List[Any]]:
    """Handle dictionary-based section specifications."""
    logger.debug("Extracting data from dict %s, %s", module, section)
    # Extract section name and column specifications
    section_name, columns = next(iter(section.items()))

    # Check if section exists
    if check_section_exists(all_module_data, section_name):
        # Find the matching section data
        section_data = next((data for key, data in all_module_data.items() if section_name in key), None)

        if section_data:
            # Convert to DataFrame and filter columns
            data = pd.DataFrame.from_dict(section_data, orient="index")
            filtered_data = filter_and_rename_columns(data, columns)
            return [filtered_data], columns

    logger.warning(f"Section '{section_name}' not found in module '{module}'")
    return [pd.DataFrame()], []
