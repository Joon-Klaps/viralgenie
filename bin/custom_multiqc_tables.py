#!/usr/bin/env python

"""Provide a command line tool to extract sequence names from cdhit's cluster files."""

import argparse
import csv
import logging
import os
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger()


class BlastConstants:
    COLUMNS = [
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
            logger.error(f"File '{fname}' doesn't end with one of {choices}")
            sys.exit(2)
        return fname_path

    parser.add_argument(
        "--clusters_summary",
        metavar="CLUSTER SUMMARY FILES",
        nargs="+",
        type=Path,
        help=" list of cluster summary files from created by the module extract clust.",
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
        "--sample_metadata",
        metavar="SAMPLE METADATA",
        help="Sample metadata file containing information on the samples, supported formats: '.csv', '.tsv'",
        type=lambda s: file_choices(("csv", "tsv"), s),
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
        logger.warning("The given input file %s is empty, it will not be used!", file)
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
        pandas.DataFrame: The concatenated dataframe.
    """
    df = pd.concat(
        [
            read_file_to_dataframe(file, **kwargs)
            for file in table_files
            if check_file_exists(file)
        ]
    )
    return df


def read_in_quast(table_files):
    """Concatenate all the cluster summary files into a single dataframe.

    Args:
        table_files (list): List of file paths to the cluster summary files.

    Returns:
        pandas.DataFrame: A dataframe containing the concatenated data from all the cluster summary files.
    """
    df = pd.DataFrame()
    if table_files:
        for file in table_files:
            if check_file_exists(file, throw_error=False):
                with open(file, "r") as f:
                    d = dict(line.strip().split("\t") for line in f)
                    df = pd.concat([df, pd.DataFrame.from_dict(d, orient="index").T])
    return df


def get_files_and_columns_of_interest(table_headers):
    """
    Get the files of interest and the columns of interest from the table headers file

    Args:
        table_headers (str): Path to the table headers file

    Returns:
        tuple: A tuple containing the files of interest (list) and the columns of interest (dic, old_colname:new_colname)
    """
    file_columns = {}
    if table_headers:
        check_file_exists(table_headers)
        # Read the yaml column annotation file for the different tables
        with open(table_headers, "r") as f:
            header_multiqc = yaml.safe_load(f)
            for bigkey, element in header_multiqc.items():
                if bigkey == "general_stats":
                    newnamekey = "multiqc"
                else:
                    newnamekey = bigkey
                columns_of_interest = {}
                if isinstance(element, list):
                    for item in element:
                        if isinstance(item, dict):
                            for key, value in item.items():
                                columns_of_interest.update(
                                    {key: f"({newnamekey}) {value}"}
                                )
                        else:
                            columns_of_interest.update(
                                {item: f"({newnamekey}) {item.replace('_',' ')}"}
                            )
                file_columns.update({bigkey: columns_of_interest})

    else:
        # Files of interest contigs:
        files_of_interest = [
            "samtools_stats",
            "umitools",
            "general_stats",
            "picard_dups",
            "ivar_variants",
            "bcftools_stats",
        ]
        file_columns = {file: {} for file in files_of_interest}

    return file_columns


def write_dataframe(df, file, comment):
    """
    Write a pandas DataFrame to a file in TSV format.

    Args:
        df (pandas.DataFrame): The DataFrame to be written.
        file (str): The file path to write the DataFrame to.
        comment (list): A list of strings to be written as comments at the beginning of the file.

    Returns:
        None
    """
    df_tsv = df.to_csv(sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC)
    with open(file, "w") as f:
        if comment:
            f.write("\n".join(comment))
            f.write("\n")
        f.write(df_tsv)


def filter_and_rename_columns(df, columns_of_interest):
    """
    Filter for columns of interest and rename those we can
    Filtering is done first on exact match and then on approximate match.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        columns_of_interest (dict): A dictionary mapping column names of interest to their desired new names.

    Returns:
        pandas.DataFrame: The filtered DataFrame with renamed columns.
    """
    if columns_of_interest:
        renamed_columns = {}
        df_columns = df.columns.tolist()

        for key, value in columns_of_interest.items():
            if key in df_columns:  # Check for an exact match first
                renamed_columns[key] = value
            else:  # If no exact match, try approximate match
                matches = [col for col in df_columns if key in col]
                if matches:
                    matched_column = matches[0]
                    renamed_columns[matched_column] = value

        # Filtering and renaming columns
        filtered_df = df[list(renamed_columns.keys())]
        filtered_df = filtered_df.rename(columns=renamed_columns, inplace=False)
        return filtered_df
    else:
        return df


def read_dataframe_from_tsv(file, **kwargs):
    """
    Read a dataframe from a tsv file.

    Args:
        file (str): The path to the file.

    Returns:
        pandas.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as table:
        df = pd.read_csv(table, sep="\t", **kwargs)
    return df


def read_dataframe_from_csv(file, **kwargs):
    """
    Read a dataframe from a csv file.

    Args:
        file (str): The path to the file.

    Returns:
        pandas.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as table:
        df = pd.read_csv(table, **kwargs)
    return df


def read_dataframe_from_yaml(file, **kwargs):
    """
    Read a dataframe from a YAML file.

    Args:
        file (str): The path to the file.

    Returns:
        pandas.DataFrame: The dataframe read from the file.
    """
    with open(file, "r") as yaml_file:
        data = yaml.safe_load(yaml_file, **kwargs)
        df = pd.DataFrame(data)
    return df


def read_file_to_dataframe(file, **kwargs):
    """
    Read a dataframe from a file.

    Args:
        file (str): The path to the file.

    Returns:
        pandas.DataFrame: The dataframe read from the file.
    """
    file_path = Path(file)
    if file_path.suffix in [".tsv", ".txt"]:  # mqc calls tsv's txts
        df = read_dataframe_from_tsv(file_path, **kwargs)
    elif file_path.suffix == ".csv":
        df = read_dataframe_from_csv(file_path, **kwargs)
    elif file_path.suffix in [".yaml", ".yml"]:
        df = read_dataframe_from_yaml(file_path, **kwargs)
    return df


def join_dataframes(df1, df2):
    """
    Join two DataFrames or a DataFrame and a Series together.

    Args:
        df1 (pd.DataFrame or pd.Series): The first DataFrame or Series to be joined.
        df2 (pd.DataFrame): The second DataFrame to be joined.

    Returns:
        pd.DataFrame: The joined DataFrame.

    """
    if isinstance(df1, pd.Series):  # Check if df1 is a Series
        df1 = df1.to_frame()  # Convert Series to DataFrame
        return pd.concat([df1, df2], axis=1, join="outer")
    elif df1.empty:
        return df2
    else:
        return df1.join(df2, how="outer")


def process_multiqc_dataframe(df):
    """
    Process the MultiQC dataframe by splitting the values in the first column by "." and setting it as the index.
    This function is required to handle the output files form quast & checkv to bring them to the same standard as MultiQC.

    Args:
        df (pandas.DataFrame): The MultiQC dataframe to be processed.

    Returns:
        pandas.DataFrame: The processed MultiQC dataframe.
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
        df (pandas.DataFrame): The dataframe containing the failed contig data.

    Returns:
        pandas.Index: The index series of the processed dataframe.
    """
    df = df.astype(str)
    df["Cluster"] = df["Cluster"].str.split(".0").str[0]
    df["Iteration"] = df["Step"].str.split(".0").str[0]
    df["id"] = df["Sample.1"] + "_" + df["Cluster"] + "_" + df["Iteration"]
    df.set_index("id", inplace=True)
    index_series = df.index
    return index_series


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
    for file_name, sub_dic in file_columns.items():
        files_of_interest = filter_files_of_interest(multiqc_data, file_name)
        if not files_of_interest:
            logger.info(
                "No files of interest were found for %s in %s", file_name, directory
            )
            continue
        df = read_dataframe_from_tsv(files_of_interest)
        df = process_dataframe(df)
        df = filter_and_rename_columns(df, sub_dic)
        multiqc_samples_df = join_dataframes(multiqc_samples_df, df)

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


def parse_annotation_data(annotation_str):
    annotation_dict = {}
    pattern = r'(?P<key>\w+)[:=]"(?P<value>[^"]+)"'
    matches = re.findall(pattern, annotation_str)
    for key, value in matches:
        annotation_dict[key] = value
    return annotation_dict


def handle_tables(table_files, header_name=False, output=False, **kwargs):
    """
    Handle multiple table files and perform concatenation and writing to output file if specified.

    Args:
        table_files (list): List of table file paths.
        header_name (bool, optional): Flag to include header name in the output file. Defaults to False.
        output (str, optional): Output file path. Defaults to False.
        **kwargs: Additional keyword arguments for concatenation.

    Returns:
        pandas.DataFrame: Concatenated table data.

    """
    result_df = pd.DataFrame()
    if table_files:
        result_df = concat_table_files(table_files, **kwargs)
    if output:
        write_dataframe(result_df, output, header_name)
    return result_df


def handle_dataframe(df, prefix, column_to_split, header=False, output=False):
    """
    Handle the given dataframe by adding a prefix to column names, splitting a specific column,
    and generating an ID based on sample, cluster, and step information.

    Args:
        df (pandas.DataFrame): The input dataframe.
        prefix (str): The prefix to add to column names.
        column_to_split (str): The column to split.
        header (bool, optional): Whether the dataframe has a header. Defaults to False.
        output (bool, optional): Whether to write the resulting dataframe to a file. Defaults to False.

    Returns:
        pandas.DataFrame: The processed dataframe.
    """
    result_df = pd.DataFrame()
    if not df.empty:
        logger.info("Handling dataframe: %s", prefix)
        df = df.add_prefix(f"({prefix}) ")

        # Apply the dynamic split function to each row in the column
        split_data = (
            df[f"({prefix}) {column_to_split}"].apply(dynamic_split).apply(pd.Series)
        )
        # take the first three columns & rename
        split_data = split_data.iloc[:, :3]
        split_data.rename(
            columns={
                split_data.columns[0]: "sample",
                split_data.columns[1]: "cluster",
                split_data.columns[2]: "step",
            },
            inplace=True,
        )

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
        df (pandas.DataFrame): The dataframe to be filtered.
        column (str): The column to filter on.
        regex_value (str): The regex value to filter on.

    Returns:
        pandas.DataFrame, pandas.DataFrame: The filtered dataframe with the regex value and the filtered dataframe without the regex value.
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
        df (pandas.DataFrame): The dataframe to drop columns from.
        columns (list): The list of columns to drop.

    Returns:
        pandas.DataFrame: The dataframe with the dropped columns.
    """
    result = df.drop(columns=[column for column in columns if column in df.columns])
    return result.copy()


def reorder_columns(df, columns):
    """
    Try to reorder columns in a dataframe and return the dataframe.

    Args:
        df (pandas.DataFrame): The dataframe to reorder columns in.
        columns (list): The list of columns to reorder.

    Returns:
        pandas.DataFrame: The dataframe with the reordered columns.
    """
    df = df[
        [column for column in columns if column in df.columns]
        + df.columns.difference(columns, sort=False).tolist()
    ]
    return df


def create_constrain_summary(df_constrain, file_columns):
    # Filter only for columns of interest
    # Some columns were already renamed, so we get the new values of them based on the original naming of mqc
    dic_columns = {
        sub_key: sub_value
        for sub_dict in file_columns.values()
        for sub_key, sub_value in sub_dict.items()
    }
    keys_to_extract = [
        "reads_mapped",
        "reads_mapped_percent",
        "number_of_SNPs",
        "10_x_pc",
        "median_coverage",
        "min_coverage",
        "max_coverage",
    ]
    columns_of_interest = [
        dic_columns[key] for key in keys_to_extract if key in dic_columns.keys()
    ]

    if not columns_of_interest:
        logger.warning(
            "No columns of interest were found to create the constrain summary table!"
        )
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
            if "segment" in df_constrain.columns
            and pd.notnull(row["species"])
            and pd.notnull(row["segment"])
            else (
                row["species"]
                if "species" in df_constrain.columns and pd.notnull(row["species"])
                else row["cluster"]
            )
        ),
        axis=1,
    )
    df_constrain = df_constrain.rename(columns={"cluster": "Constrain id"})

    # Remove columns that are not needed anymore
    df_constrain = drop_columns(df_constrain, ["species", "segment"])

    # Convert dataframe to long and then extra wide
    df_long = df_constrain.melt(
        id_vars=["idgroup", "sample name"], var_name="variable", value_name="Value"
    )
    # Remove rows with NaN values & duplicates
    df_long = df_long.dropna()
    df_long = df_long.drop_duplicates()
    df_long["grouped variable"] = df_long["idgroup"] + " - " + df_long["variable"]
    df_long.drop(columns=["idgroup", "variable"], inplace=True)
    # Convert to wide format
    df_wide = df_long.pivot(
        index=["sample name"], columns="grouped variable", values="Value"
    )
    df_wide.reset_index(inplace=True)

    return df_wide


def main(argv=None):
    """
    Main function for creating custom tables for MultiQC.

    Args:
        argv (list): List of command line arguments.

    Returns:
        int: Exit code.
    """
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # General stats - Cluster summariesx
    if args.clusters_summary:
        cluster_header = get_header(args.comment_dir, "clusters_summary_mqc.txt")
        handle_tables(args.clusters_summary, cluster_header, "summary_clusters_mqc.tsv")

    # General Stats - Sample metadata
    if args.sample_metadata:
        sample_header = get_header(args.comment_dir, "sample_metadata_mqc.txt")
        handle_tables([args.sample_metadata], sample_header, "sample_metadata_mqc.tsv")

    # CLuster table - Checkv summary
    checkv_df = handle_tables(args.checkv_files)
    checkv_header = []
    if args.save_intermediate:
        checkv_header = get_header(args.comment_dir, "checkv_mqc.txt")
    if not checkv_df.empty:
        checkv_df = handle_dataframe(
            checkv_df, "checkv", "contig_id", checkv_header, "summary_checkv_mqc.tsv"
        )

    # CLuster table - Quast summary
    quast_df = read_in_quast(args.quast_files)
    quast_header = []
    if args.save_intermediate:
        quast_header = get_header(args.comment_dir, "quast_mqc.txt")
    if not quast_df.empty:
        quast_df = handle_dataframe(
            quast_df, "quast", "Assembly", quast_header, "summary_quast_mqc.tsv"
        )
        # Most of the columns are not good for a single contig evaluation
        quast_df["(quast) # N's"] = (
            pd.to_numeric(quast_df["(quast) # N's per 100 kbp"])
            * pd.to_numeric(quast_df["(quast) Largest contig"])
            / 100000
        )
        quast_df = quast_df.astype({"(quast) # N's": int})
        quast_df["(quast) % N's"] = round(
            pd.to_numeric(quast_df["(quast) # N's per 100 kbp"]) / 1000, 2
        )
        quast_df = quast_df[
            ["(quast) # N's", "(quast) % N's", "(quast) # N's per 100 kbp"]
        ]

    # CLuster table - Blast summary
    blast_df = handle_tables(args.blast_files, header=None)
    blast_header = []
    if args.save_intermediate:
        blast_header = get_header(args.comment_dir, "blast_mqc.txt")
    if not blast_df.empty:
        # Read the blast summary file
        blast_df.columns = BlastConstants.COLUMNS

        # Filter on best hit per contig and keep only the best hit
        blast_df = blast_df.sort_values("bitscore", ascending=False).drop_duplicates(
            "query"
        )

        blast_df = handle_dataframe(
            blast_df, "blast", "query", blast_header, "summary_blast_mqc.tsv"
        )

    # CLuster table - mmseqs easysearch summary - annotation section
    annotation_df = handle_tables(args.annotation_files, header=None)
    if not annotation_df.empty:
        annotation_df.columns = BlastConstants.COLUMNS
        # Filter on best hit per contig and keep only the best hit
        annotation_df = annotation_df.sort_values(
            "bitscore", ascending=False
        ).drop_duplicates("query")

        # Extract all key-value pairs into separate columns
        df_extracted = (
            annotation_df["subject title"].apply(parse_annotation_data).apply(pd.Series)
        )

        # Concatenate the original DataFrame with the extracted columns
        annotation_df = pd.concat([annotation_df, df_extracted], axis=1)

        # Remove the blast columns (but not query), we only want annotation data (but not genome_name)
        annotation_df = drop_columns(annotation_df, BlastConstants.COLUMNS[1:])
        annotation_df = handle_dataframe(
            annotation_df, "annotation", "query", blast_header, "summary_anno_mqc.tsv"
        )
        annotation_df["(annotation) taxon_id"] = annotation_df[
            "(annotation) taxon_id"
        ].astype(str)

    # CLuster table - Blast summary
    annotation_df = handle_tables(args.annotation_files, header=None)
    if not annotation_df.empty:
        annotation_df.columns = BlastConstants.COLUMNS
        # Filter on best hit per contig and keep only the best hit
        annotation_df = annotation_df.sort_values(
            "bitscore", ascending=False
        ).drop_duplicates("query")

        # Extract all key-value pairs into separate columns
        df_extracted = (
            annotation_df["subject title"].apply(parse_annotation_data).apply(pd.Series)
        )

        # Concatenate the original DataFrame with the extracted columns
        annotation_df = pd.concat([annotation_df, df_extracted], axis=1)

        # Remove the blast columns (but not query), we only want annotation data (but not genome_name)
        annotation_df = drop_columns(annotation_df, BlastConstants.COLUMNS[1:])
        annotation_df = handle_dataframe(
            annotation_df, "annotation", "query", blast_header, "summary_anno_mqc.tsv"
        )
        annotation_df["(annotation) taxon_id"] = annotation_df[
            "(annotation) taxon_id"
        ].astype(str)

    # CLuster table -  Multiqc output txt files
    if args.multiqc_dir:
        # Check if the given files exist
        check_file_exists(args.multiqc_dir)

        # Extract files & columns specified by the user through the table headers file
        file_columns = get_files_and_columns_of_interest(args.table_headers)

        # Read the multiqc data txt files & select & rename
        multiqc_contigs_df = read_data(
            args.multiqc_dir, file_columns, process_multiqc_dataframe
        )

        # If we are empty, just quit
        if multiqc_contigs_df.empty:
            logger.warning(
                "No data was found from MULTIQC to create the contig overview table!"
            )
            return 0

        # Write the complete dataframe to a file
        logger.info("Writing intermediate file: contigs_intermediate.tsv")
        write_dataframe(
            multiqc_contigs_df.reset_index(inplace=False),
            "contigs_intermediate.tsv",
            [],
        )

        # Join with the custom contig tables
        logger.info("Joining dataframes")
        multiqc_contigs_df = multiqc_contigs_df.join(checkv_df, how="outer")
        multiqc_contigs_df = multiqc_contigs_df.join(quast_df, how="outer")
        multiqc_contigs_df = multiqc_contigs_df.join(blast_df, how="outer")
        multiqc_contigs_df = multiqc_contigs_df.join(annotation_df, how="outer")

        # adding a tag saying that contig faild qc check
        logger.info("Adding failed contig QC check")
        failed_contigs = {"failed_mapped": {}, "failed_contig_quality": {}}
        failed_contig_df = read_data(
            args.multiqc_dir, failed_contigs, process_failed_contig_dataframe
        )
        if not failed_contig_df.empty:
            multiqc_contigs_df["Contig failed QC check"] = (
                multiqc_contigs_df.index.isin(failed_contig_df)
            )

        # If we are empty, just quit
        if multiqc_contigs_df.empty:
            logger.warning("No data was found to create the contig overview table!")
            return 0

        # Keep only those rows we can split up in sample, cluster, step
        logger.info("Splitting up the index column in sample name, cluster, step")
        mqc_contigs_sel = multiqc_contigs_df.reset_index().rename(
            columns={multiqc_contigs_df.index.name: "index"}
        )
        mqc_contigs_sel = mqc_contigs_sel[
            mqc_contigs_sel["index"].str.contains("_", na=False)
        ]
        mqc_contigs_sel[["sample name", "cluster", "step"]] = mqc_contigs_sel[
            "index"
        ].str.split("_", n=3, expand=True)

        # Reorder the columns
        logger.info("Reordering columns")
        final_columns = (
            ["index", "sample name", "cluster", "step"]
            + [column for column in mqc_contigs_sel.columns if "annotation" in column]
            + [column for column in mqc_contigs_sel.columns if "blast" in column]
            + [column for column in mqc_contigs_sel.columns if "checkv" in column]
            + [column for column in mqc_contigs_sel.columns if "QC check" in column]
            + [column for column in mqc_contigs_sel.columns if "quast" in column]
        )
        mqc_contigs_sel = reorder_columns(mqc_contigs_sel, final_columns)

        # split up denovo constructs and mapping (-CONSTRAIN) results
        logger.info("Splitting up denovo constructs and mapping (-CONSTRAIN) results")
        contigs_mqc, constrains_mqc = filter_constrain(
            mqc_contigs_sel, "cluster", "-CONSTRAIN"
        )

        # Write the final dataframe to a file
        logger.info("Writing Denovo constructs table file: contigs_overview_mqc.tsv")
        header_clusters_overview = get_header(
            args.comment_dir, "contig_overview_mqc.txt"
        )
        write_dataframe(
            contigs_mqc, "contigs_overview_mqc.tsv", header_clusters_overview
        )

        # Separate table for mapping constrains
        if not constrains_mqc.empty:
            header_mapping_seq = get_header(
                args.comment_dir, "mapping_constrains_mqc.txt"
            )

            # Add constrain metadata to the mapping constrain table
            constrain_meta = handle_tables([args.mapping_constrains])
            # drop unwanted columns & reorder
            constrain_meta = drop_columns(constrain_meta, ["sequence", "samples"])
            constrains_mqc = constrains_mqc.merge(
                constrain_meta, how="left", left_on="cluster", right_on="id"
            )
            constrains_mqc = reorder_columns(
                constrains_mqc,
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
            logger.info("Writing mapping long table: mapping_constrains_mqc.tsv")
            write_dataframe(
                constrains_mqc, "mapping_constrains_mqc.tsv", header_mapping_seq
            )

            # add mapping summary to sample overview table in ... wide format with species & segment combination
            logger.info("Creating mapping constrain summary (wide) table")
            constrains_summary_mqc = create_constrain_summary(
                constrains_mqc, file_columns
            )
            if not constrains_summary_mqc.empty:
                header_mapping_summary = get_header(
                    args.comment_dir, "mapping_constrains_summary_mqc.txt"
                )
                write_dataframe(
                    constrains_summary_mqc,
                    "mapping_constrains_summary_mqc.tsv",
                    header_mapping_summary,
                )
    return 0


if __name__ == "__main__":
    sys.exit(main())
