#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.


"""Provide a command line tool to create several custom mqc report files."""

from __future__ import annotations

import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from utils.constant_variables import FILES_OF_INTEREST

# Optional dependencies - set to None if not available
try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import yaml
except ImportError:
    yaml = None

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

logger = logging.getLogger()


def _check_pandas_available():
    """Check if pandas is available, raise ImportError if not."""
    if pd is None:
        raise ImportError(
            "pandas is required for this function but is not installed. "
            "Please install it with: pip install pandas"
        )


def _check_yaml_available():
    """Check if PyYAML is available, raise ImportError if not."""
    if yaml is None:
        raise ImportError(
            "PyYAML is required for this function but is not installed. "
            "Please install it with: pip install pyyaml"
        )


def _check_biopython_available():
    """Check if Biopython is available, raise ImportError if not."""
    if SeqIO is None:
        raise ImportError(
            "Biopython is required for this function but is not installed. "
            "Please install it with: pip install biopython"
        )


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
    _check_pandas_available()
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
    _check_pandas_available()
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
    _check_pandas_available()

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
    _check_pandas_available()

    file_path = Path(file)
    if os.path.getsize(file_path) == 0:
        logger.debug("File is empty %s", file_path)
        return pd.DataFrame()

    # pandas decompresses on the fly, we only need to look past the .gz to find the format
    suffixes = file_path.suffixes
    suffix = suffixes[-2] if suffixes and suffixes[-1] == ".gz" and len(suffixes) > 1 else file_path.suffix

    if suffix in [
        ".tsv",
        ".txt",
    ]:  # mqc calls tsv's txts, bed files are gzipped
        df = df_from_tsv(file_path, **kwargs)
    elif suffix == ".csv":
        df = df_from_csv(file_path, **kwargs)
    elif suffix in [".yaml", ".yml"]:
        df = df_from_yaml(file_path, **kwargs)
    elif suffix in [".json"]:
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
    _check_pandas_available()

    df = pd.read_csv(file, sep="\t", **kwargs)
    return df


def df_from_csv(file: str, **kwargs) -> pd.DataFrame:
    """
    Read a dataframe from a csv file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    _check_pandas_available()

    df = pd.read_csv(file, **kwargs)
    return df


def df_from_yaml(file: str, **kwargs) -> pd.DataFrame:
    import yaml
    """
    Read a dataframe from a YAML file.

    Args:
        file (str): The path to the file.

    Returns:
        pd.DataFrame: The dataframe read from the file.
    """
    _check_yaml_available()
    _check_pandas_available()

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
    _check_pandas_available()

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
    _check_pandas_available()

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

    _check_yaml_available()
    check_file_exists(table_headers)
    yaml_data = yaml.safe_load(table_headers.read_text())

    return yaml_data


def index_fasta(fasta_path, format="fasta"):
    """
    Create an indexed dictionary of sequences from a FASTA file.

    Tries SeqIO.index() first (memory-efficient), falls back to SeqIO.to_dict()
    with duplicate removal if duplicate keys are found.

    Args:
        fasta_path: Path to the FASTA file
        format: Sequence format (default: "fasta")

    Returns:
        Dict-like object mapping sequence IDs to SeqRecord objects

    Raises:
        ValueError: If duplicates found AND file is too large to fit in memory
    """
    _check_biopython_available()

    try:
        return SeqIO.index(str(fasta_path), format)
    except ValueError as e:
        if "Duplicate key" not in str(e):
            raise

        logger.warning(
            "Duplicate sequence IDs found in %s: %s. "
            "Falling back to in-memory indexing (keeping first occurrence).",
            fasta_path,
            e,
        )

        try:
            return _to_dict_first_occurrence(SeqIO.parse(str(fasta_path), format))
        except MemoryError:
            raise ValueError(
                f"Duplicate sequence IDs found in {fasta_path} and file is too large "
                f"to load into memory. Please deduplicate the FASTA file first "
                f"(e.g., using 'seqkit rmdup -n'). Original error: {e}"
            ) from e


def _to_dict_first_occurrence(sequences) -> dict:
    """
    Convert sequences to dict, keeping only FIRST occurrence of duplicate IDs.

    Args:
        sequences: Iterable of SeqRecord objects

    Returns:
        Dict mapping sequence IDs to SeqRecord objects
    """
    result = {}
    duplicates = []
    for record in sequences:
        if record.id not in result:
            result[record.id] = record
        else:
            duplicates.append(record.id)

    if duplicates:
        unique_dups = list(dict.fromkeys(duplicates))  # Preserve order, remove dups
        logger.warning(
            "Skipped %d duplicate sequence(s), keeping first occurrence: %s%s",
            len(duplicates),
            ", ".join(unique_dups[:5]),
            "..." if len(unique_dups) > 5 else "",
        )

    return result
