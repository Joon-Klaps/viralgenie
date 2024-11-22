#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import logging
from typing import Dict, List, Union

import pandas as pd

logger = logging.getLogger()


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


def fill_group_na(s):
    return s.infer_objects(copy=False).ffill().bfill()


def dynamic_split(index_str):
    """
    Split the index string into sample name, cluster, and step.
    """
    parts = index_str.split("_")
    return parts[0], parts[1], "_".join(parts[2:])


def join_df(base_df: pd.DataFrame, dfs: List[pd.DataFrame]) -> pd.DataFrame:
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
