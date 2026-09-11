#!/usr/bin/env python

# Originally written by Joon Klaps and released under the MIT license.
# See git repository (https://github.com/nf-core/viralmetagenome) for full license text.


import logging
import re
from typing import Dict, List, Union, Tuple, Optional, Any

import pandas as pd

from utils.constant_variables import BLAST_COLUMNS, CONSTRAINT_GENERAL_STATS_COLUMNS, COLUMN_MAPPING, READ_DECLARATION
from utils.file_tools import filelist_to_df, read_file_to_df
from utils.pandas_tools import (
    coalesce_constraint,
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


def process_annotation_df(annotation_df, metadata_df=None):
    """
    Process the annotation DataFrame.

    Args:
        annotation_df (pd.DataFrame): The annotation DataFrame.
        metadata_df (pd.DataFrame): Optional annotation metadata table, as returned by
            read_annotation_metadata(). When given, the annotation fields are looked up
            in it instead of being parsed out of the fasta header.

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
        annotation_df = extract_annotation_data(annotation_df, metadata_df)

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


def read_annotation_metadata(metadata_file) -> pd.DataFrame:
    """
    Read an annotation metadata table that describes the sequences of the annotation database.

    The first column must hold the sequence identifiers as they appear in the fasta header
    (the part of the '>' line up to the first whitespace); every other column becomes an
    annotation field. Supported formats are csv, tsv and their gzipped counterparts.

    Args:
        metadata_file: Path to the metadata table, or None.

    Returns:
        pd.DataFrame: The metadata indexed by sequence identifier, empty if no file was given.
    """
    if not metadata_file:
        return pd.DataFrame()

    metadata_df = read_file_to_df(metadata_file, dtype=str)
    if metadata_df.empty:
        logger.warning("The annotation metadata file %s is empty, ignoring it.", metadata_file)
        return metadata_df

    if len(metadata_df.columns) < 2:
        logger.warning(
            "The annotation metadata file %s has no annotation columns next to the identifier column, ignoring it.",
            metadata_file,
        )
        return pd.DataFrame()

    id_column = metadata_df.columns[0]
    logger.info(
        "Reading annotation metadata from %s: using '%s' as the sequence identifier and %s as annotation fields.",
        metadata_file,
        id_column,
        ", ".join(f"'{column}'" for column in metadata_df.columns[1:]),
    )

    metadata_df = metadata_df.set_index(id_column)
    duplicates = metadata_df.index.duplicated()
    if duplicates.any():
        logger.warning(
            "The annotation metadata file %s has %d duplicated identifiers, keeping the first occurrence of each.",
            metadata_file,
            int(duplicates.sum()),
        )
        metadata_df = metadata_df[~duplicates]

    return metadata_df


def extract_annotation_data(df, metadata_df=None):
    """
    Add the annotation fields of the best hit as separate columns.

    The fields come from the metadata table when one is given, and are parsed out of the
    fasta header of the hit otherwise. An identifier is matched as it is written first and
    on its versionless form after that ('NC_078521.1' -> 'NC_078521'), so that a database
    and a metadata table that disagree on the suffix still line up.
    """
    if metadata_df is None or metadata_df.empty:
        extracted = df["subject title"].astype(str).apply(parse_annotation_data).apply(pd.Series)
    else:
        # exact identifiers come first so that they win over the versionless ones
        metadata_df = metadata_df.rename(index=str)
        lookup = pd.concat([metadata_df, metadata_df.rename(index=lambda i: i.split(".")[0])])
        lookup = lookup[~lookup.index.duplicated()]
        # identifiers that look numeric ('754189.6') are read as floats, the lookup needs the text
        subject = df["subject"].astype(str)
        keys = subject.where(subject.isin(lookup.index), subject.str.split(".").str[0])
        extracted = lookup.reindex(keys)
        extracted.index = df.index
        # counted on the keys, a row that was found can have empty annotation values of its own
        missing = int((~keys.isin(lookup.index)).sum())
        if missing:
            logger.warning(
                "%d of the %d annotated contigs hit a sequence that is not in the annotation metadata file.",
                missing,
                len(df),
            )

    # An annotation field that is named like one of the search result columns would end up
    # as a duplicated column and break every lookup on it further down. Names are compared
    # without case, as columns differing only in case ('length' from the search and 'Length'
    # from a metadata table) are read as one by plenty of table readers.
    clashes = extracted.columns[extracted.columns.str.casefold().isin(df.columns.str.casefold())]
    if not clashes.empty:
        logger.warning(
            "Ignoring annotation field(s) %s as they collide with the search result columns.",
            ", ".join(clashes),
        )

    return pd.concat([df, extracted.drop(columns=clashes)], axis=1)


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
    df = df.loc[:, ~df.columns.duplicated()]
    return reorder_columns(df.dropna(subset=["step"]), list(dict.fromkeys(final_columns)))


def filter_constraint(dataframe, column, value):
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
    locations = df[column].str.contains(value) | df["step"].str.contains("constraint")

    # Filter
    df_with_value = df[locations]
    df_without_value = df[~locations]
    # Remove from column
    df_with_value.loc[:, column] = df_with_value[column].str.replace(value, "")
    df_with_value.loc[:, "index"] = df_with_value["index"].str.replace(value, "")

    return df_without_value.dropna(axis=1, how="all"), df_with_value.dropna(axis=1, how="all")

def get_anchor(mqc: object, module: str) -> str | None:
    for m in mqc.report.modules:
        if m.name.lower() == module.lower():
            return m.anchor
    return None

def order_columns_by_idgroup(df_wide: pd.DataFrame, original_columns: List[str]) -> List[str]:
    """
    Order columns by species-segment group (idgroup) while maintaining original column order within each group.

    Args:
        df_wide (pd.DataFrame): The wide-format DataFrame with columns formatted as "idgroup - column"
        original_columns (List[str]): The original column order to maintain within each idgroup

    Returns:
        List[str]: Ordered list of column names
    """
    # Get all unique species-segment combinations (idgroups)
    idgroups = sorted(set([col.split(' - ')[0] for col in df_wide.columns]))

    # Create an ordered column list that groups by species-segment while maintaining original column order for each group
    ordered_columns = []
    for idgroup in idgroups:
        # For each species-segment group, add columns in the original order
        for col in original_columns:
            group_col = f"{idgroup} - {col}"
            if group_col in df_wide.columns:
                ordered_columns.append(group_col)

    return ordered_columns

def create_constraint_summary(df_constraint: pd.DataFrame, file_columns: List[Union[str, Dict[str, str]]]) -> pd.DataFrame:
    """
    Create a summary table for the constraint data.

    Args:
        df_constraint (pd.DataFrame): The constraint data DataFrame.
        file_columns (List): A columns dictionary with old names & new names.

    Returns:
        pd.DataFrame: The constraint summary table.
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

    columns_of_interest = [dic_columns.get(key, key) for key in CONSTRAINT_GENERAL_STATS_COLUMNS]

    logger.debug("columns_of_interest: %s", columns_of_interest)

    logger.debug("available columns: %s", df_constraint.columns)

    if not columns_of_interest:
        logger.warning("No columns of interest were found to create the constraint summary table!")
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

    df_columns = df_constraint.columns.tolist()

    present_columns = []
    for name in columns_of_interest:
        if name in df_columns:  # Check for an exact match first
            present_columns.append(name)
        else:  # If no exact match, try approximate match
            matches = [col for col in df_columns if name in col]
            if matches:
                matched_column = matches[0]
                present_columns.append(matched_column)

    df_constraint = df_constraint[present_columns]

    if df_constraint.empty:
        logger.warning("The constraint DataFrame is empty.")
        return df_constraint

    df_constraint = df_constraint.rename(columns=COLUMN_MAPPING)

    # Reformat dataframe to long based on following:
    #   Species & Segment
    #   Species
    #   ID (Cluster)
    df_constraint.loc[:, "idgroup"] = df_constraint.apply(
        lambda row: (
            f"{row['species']} ({row['segment']})"
            if "segment" in df_constraint.columns and pd.notnull(row["species"]) and pd.notnull(row["segment"])
            else (row["species"] if "species" in df_constraint.columns and pd.notnull(row["species"]) else row["cluster"])
        ),
        axis=1,
    )
    df_constraint = df_constraint.rename(columns={"cluster": "Constraint id"})

    # Remove columns that are not needed anymore
    df_constraint = drop_columns(df_constraint, ["species", "segment"])

    # Get original column order from df_constraint (excluding sample)
    original_columns = df_constraint.columns.drop(['sample']).tolist()

    # Convert dataframe to long and then extra wide
    df_long = df_constraint.melt(id_vars=["idgroup", "sample"], var_name="variable", value_name="Value")
    # Remove rows with NaN values & duplicates
    df_long = df_long.dropna()
    df_long = df_long.drop_duplicates()
    df_long["grouped variable"] = df_long["idgroup"] + " - " + df_long["variable"]
    df_long.drop(columns=["idgroup", "variable"], inplace=True)

    # Convert to wide format
    df_wide = df_long.pivot(index=["sample"], columns="grouped variable", values="Value")

    # Reorder the columns by species-segment groups while maintaining original column order
    ordered_columns = order_columns_by_idgroup(df_wide, original_columns)
    df_wide = df_wide[ordered_columns]

    df_wide.reset_index(inplace=True)

    return df_wide


def reformat_constraint_df(df, columns, args):
    """
    Reformat the constraint dataframe.
    """
    # Separate table for mapping constraints
    if df.empty:
        logger.warning("The constraint DataFrame is empty.")
        return df, df

    # Add constraint metadata to the mapping constraint table
    constraint_meta = filelist_to_df([args.mapping_constraints])

    # drop unwanted columns & reorder
    constraint_meta = drop_columns(constraint_meta, ["sequence", "samples"])
    df = df.merge(constraint_meta, how="left", left_on="cluster", right_on="id")
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
    logger.info("Creating mapping constraint summary (wide) table")
    mapping_constraints_summary = create_constraint_summary(df, columns).set_index("sample")

    logger.info("Coalescing columns")
    coalesced_constraints = coalesce_constraint(df)
    coalesced_constraints = drop_columns(coalesced_constraints, ["id", "selection", "rank"])
    return coalesced_constraints, mapping_constraints_summary


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


def extract_mqc_from_str_section(all_module_data: Dict, section: Optional[str], module: str) -> Tuple[List[pd.DataFrame], List[Any]]:
    """Handle simple string or None section cases."""
    logger.debug("Extracting data from simple str %s", module)
    if not section:
        # Return all data if no specific section is specified
        return [pd.DataFrame.from_dict(all_module_data, orient="index")], []

    # Check if specific section exists
    if check_section_exists(all_module_data, section):
        logger.info("Str-Section: %s found of module %s", section, module)
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
            subsection_dfs, subsection_columns = extract_mqc_from_str_section(all_module_data, subsection, module)
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
        logger.info("Dict-Section: %s found of module %s", section_name, module)
        # Find the matching section data
        section_data = next((data for key, data in all_module_data.items() if section_name in key), None)

        if section_data:
            # Convert to DataFrame and filter columns
            data = pd.DataFrame.from_dict(section_data, orient="index")
            filtered_data = filter_and_rename_columns(data, columns)
            return [filtered_data], columns

    logger.warning(f"Section '{section_name}' not found in module '{module}'")
    return [pd.DataFrame()], []

def get_read_suffix(namespace: str, title: str) -> str | None:
    """
    Get the read suffix based on the namespace and title.

    Args:
        namespace: The namespace to check
        title: The title to lookup

    Returns:
        str | None: The read suffix if found, None otherwise
    """
    title = title.strip()
    if title not in READ_DECLARATION:
        return None

    title_config = READ_DECLARATION[title]
    if any(pattern in namespace for pattern in title_config['namespace_patterns']):
        return title_config['suffix']

    return None
