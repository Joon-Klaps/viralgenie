#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import multiqc as mqc
import numpy as np
import pandas as pd
from multiqc.plots import bargraph, table
from multiqc.types import Anchor
from utils.constant_variables import CLUSTER_PCONFIG
from utils.file_tools import filelist_to_df, get_module_selection, read_in_quast, write_df
from utils.module_data_processing import (
    filter_constrain,
    generate_ignore_samples,
    generate_indexed_df,
    process_annotation_df,
    process_blast_df,
    reformat_constrain_df,
    reformat_custom_df,
)
from utils.pandas_tools import (
    filter_and_rename_columns,
    join_df,
    select_columns,
)

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
        blast_df = process_blast_df(blast_df)

    # Cluster table - mmseqs easysearch summary (annotation section)
    annotation_df = filelist_to_df(args.annotation_files, header=None)
    if not annotation_df.empty:
        annotation_df = process_annotation_df(annotation_df)

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


def handle_module_data(
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
                    subsection_data, columns = handle_module_data(module, subsection)
                    if isinstance(subsection_data, pd.DataFrame):
                        result_df.extend(subsection_data)
                    if isinstance(columns, list):
                        result_list.extend(columns)
                return result_df, result_list
    elif isinstance(section, dict):
        # We just have {section: List[Union[column_name, Dict[column_name, column_rename]]}
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
            module_data, columns = handle_module_data(module, section)
            logger.debug("Data for %s: %s", module, module_data)

        data.extend(module_data)
        columns_result.extend(columns)

    logger.debug("Data list: %s", data)

    return join_df(result, data) if data else result, columns_result


def write_results(contigs_mqc, constrains_mqc, constrains_genstats, args) -> int:
    """
    Write the results to files.
    """

    if not contigs_mqc.empty:
        logger.info("Writing Unfiltered Denovo constructs table file: contigs_all.tsv")
        write_df(contigs_mqc, "contigs_all.tsv", [])
        contigs_mqc.set_index("index", inplace=True)
        table_plot = contigs_mqc[~contigs_mqc.index.isin(generate_ignore_samples(contigs_mqc))]
        mqc.add_custom_content_section(
            name="Denovo Construct Overview",
            anchor=Anchor("contigs_all"),
            description="The table below shows the overview of the denovo constructs with refinement.",
            plot=table.plot(data=table_plot.to_dict(orient="index")),
        )

    if not constrains_mqc.empty:
        logger.info("Writing Unfiltered Mapping constructs table file: mapping_all.tsv")
        write_df(constrains_mqc, "mapping_all.tsv", [])
        constrains_mqc.set_index("index", inplace=True)
        mqc.add_custom_content_section(
            name="Mapping Construct Overview",
            anchor=Anchor("mapping_all"),
            description="The table below shows the overview of the mapping constructs with refinement.",
            plot=table.plot(data=constrains_mqc.to_dict(orient="index")),
        )

    if not constrains_genstats.empty:
        write_df(constrains_genstats, "mapping_constrains_summary.tsv", [])
        # Add to mqc
        module = mqc.BaseMultiqcModule(name="Mapping Constrains Summary", anchor=Anchor("custom_data"))
        content = constrains_genstats.to_dict(orient="index")
        module.general_stats_addcols(content)
        mqc.report.modules.append(module)

    # TODO correctly insert metadata of:
    #   -  versions
    #   -  Not all mapping data is in the general stats table, while it should be
    #   -  Double check for any other loss of information.
    mqc.write_report(
        make_data_dir=True,
        data_format="tsv",
        export_plots=False,
        force=True,
         )

    return 0


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
    mqc_custom_df = join_df(mqc_custom_df, custom_tables)

    if mqc_custom_df.empty:
        logger.warning("No data was found to create the contig overview table!")
        return 0

    # 3.3 reformat the dataframe
    mqc_custom_df = reformat_custom_df(mqc_custom_df)

    # 3.4 split up denovo constructs and mapping (-CONSTRAIN) results
    logger.info("Splitting up denovo constructs and mapping (-CONSTRAIN) results")
    contigs_mqc, constrains_mqc = filter_constrain(mqc_custom_df, "cluster", "-CONSTRAIN")

    coalesced_constrains, constrains_genstats = reformat_constrain_df(constrains_mqc, renamed_columns, args)

    write_results(contigs_mqc, coalesced_constrains, constrains_genstats, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
