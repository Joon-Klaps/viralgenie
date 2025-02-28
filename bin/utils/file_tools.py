#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import yaml

from utils.constant_variables import FILES_OF_INTEREST

logger = logging.getLogger()


def check_file_exists(file: str, throw_error=True) -> bool:
    """Check if the given files exist.

    Args:
        file (str): The path to the file to be checked.
        throw_error (bool, optional): Whether to throw an error and exit the program if the file is not found.
            Defaults to True.

    Returns:
        bool: True if the file exists and is not empty, False otherwise.
    """
    if not file:
        logger.warning("The given input file %s was not a file!", file)
        return False
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


def concat_table_files(table_files: List[str], **kwargs) -> pd.DataFrame:
    """Concatenate all the cluster summary files into a single dataframe.

    Args:
        table_files (list): List of file paths to be concatenated.
        **kwargs: Additional keyword arguments to be passed to pd.read_csv().

    Returns:
        pd.DataFrame: The concatenated dataframe.
    """
    try:
        valid_dfs = [read_file_to_df(file, **kwargs) for file in table_files if check_file_exists(file)]

        if not valid_dfs:
            logging.warning("Warning concatenating files: %s", table_files)
            logging.warning("No valid files found to concatenate.")
            return pd.DataFrame()

        df = pd.concat(valid_dfs)
        return df

    except ValueError as e:
        logging.warning("Error concatenating files: %s\n%s", table_files, e)
        return pd.DataFrame()


def read_in_quast(table_files: List[str]) -> pd.DataFrame:
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


def write_df(df: pd.DataFrame, file: str, comment: Optional[List[str]] = None) -> None:
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


def read_file_to_df(file: str, **kwargs) -> pd.DataFrame:
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


def df_from_tsv(file: str, **kwargs) -> pd.DataFrame:
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


def df_from_csv(file: str, **kwargs) -> pd.DataFrame:
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


def df_from_yaml(file: str, **kwargs) -> pd.DataFrame:
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


def df_from_json(file: str, **kwargs) -> pd.DataFrame:
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
            data["filename"] = filename + "_constraint"
        df = pd.DataFrame([data])
    return df


def filelist_to_df(table_files: List[str], header_name: Optional[List[str]] = None, output: Optional[str] = None, **kwargs) -> pd.DataFrame:
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
        write_df(result_df, output, header_name)
    return result_df


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
