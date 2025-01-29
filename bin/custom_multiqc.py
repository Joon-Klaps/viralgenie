#!/usr/bin/env python

"""Provide a command line tool to create several custom mqc report files."""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import multiqc as mqc
import pandas as pd
from multiqc.plots import bargraph
from multiqc.types import Anchor
from utils.constant_variables import CLUSTER_PCONFIG
from utils.file_tools import filelist_to_df, get_module_selection, read_in_quast, write_df
from utils.module_data_processing import *
from utils.pandas_tools import filter_and_rename_columns, join_df, reorder_columns, select_columns

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
        "--mapping_constraints",
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

def get_failed_samples(samples: List[str]) -> List[str]:
    """
    Get failed samples from the modules
        - sample_low_reads
        - samples_without_contigs
    """
    if (samples_low_reads :=  get_module_data(mqc, 'samples_low_reads')):
        logger.info("samples_low_reads %s", samples_low_reads)
        samples.extend([k for k in samples_low_reads.keys()])

    if (samples_without_contigs := get_module_data(mqc, 'samples_without_contigs')):
        logger.info("samples_without_contigs %s", samples_without_contigs)
        samples.extend([k for k in samples_without_contigs.keys() ])

    return samples

def load_custom_data(args) -> List[pd.DataFrame]:
    """
    Load custom data from files and process it to a list of dataframes.
    """
    result = []
    # Clusters overview - mini multiqc module
    clusters_summary_df = filelist_to_df(args.clusters_summary)
    if not clusters_summary_df.empty:
        clusters_summary_df.set_index("Sample name", inplace=True)
        # Adding to general stats
        module = mqc.BaseMultiqcModule(name="Cluster Summary", anchor=Anchor("cluster-summary"))
        module.general_stats_addcols(clusters_summary_df.to_dict(orient="index"))
        # Custom barplot -  Clusters sample
        plot_df = select_columns(clusters_summary_df, ["# Clusters", "# Removed clusters"])
        plot = bargraph.plot(data=plot_df.to_dict(orient="index"), pconfig=CLUSTER_PCONFIG)
        module.add_section(
            anchor=Anchor("cluster-summary"), plot=plot, description="Number of identified contig clusters per sample after assembly."
        )
        mqc.report.modules.append(module)

    # General Stats - Sample metadata
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
        result.extend([checkv_df])

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
        result.extend([quast_df])

    # Cluster table - Blast summary
    blast_df = filelist_to_df(args.blast_files, header=None)
    if not blast_df.empty:
        blast_df = process_blast_df(blast_df)
        result.extend([blast_df])

    # MASH screen used for reference selection summarisation
    screen_df = filelist_to_df(args.screen_files)
    if not screen_df.empty:
        screen_df = generate_indexed_df(screen_df, "mash-screen", "filename")
        screen_df = screen_df.astype(str)
        result.extend([screen_df])

    # Cluster table - mmseqs easysearch summary (annotation section)
    annotation_df = filelist_to_df(args.annotation_files, header=None)
    if not annotation_df.empty:
        annotation_df = process_annotation_df(annotation_df)
        result.extend([annotation_df])

    # Cluster table - cluster summary of members & centroids
    clusters_df = filelist_to_df(args.clusters_files)
    if not clusters_df.empty:
        clusters_df = clusters_df.add_prefix("(cluster) ")
        clusters_df = clusters_df.rename(columns={"(cluster) sample": "sample", "(cluster) cluster": "cluster"})

    return result, clusters_df


def get_general_stats_data_mod(sample: Optional[str] = None) -> Dict:
    """
    Return parsed general stats data, indexed by sample, then by data key. If sample is specified,
    return only data for that sample.

    @param sample: Sample name
    @return: Dict of general stats data indexed by sample and data key
    """

    data: Dict[str, Dict] = defaultdict(dict)
    for rows_by_group, header in zip(mqc.report.general_stats_data, mqc.report.general_stats_headers):
        for s, rows in rows_by_group.items():
            if sample and s != sample:
                continue
            for row in rows:
                for key, val in row.data.items():
                    if key in header:
                        namespace = header[key].get("namespace", key).replace("SAMPLE: ", "")
                        final_key = f"{namespace}. {header[key].get('title', key)}" if header[key].get("title") else key
                        data[s][final_key] = val
    if sample:
        if not data:
            return {}
        return data[sample]

    return data


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


def extract_module_data(
    module: str, section: Union[str, None, List[Union[str, Dict[str, str]]], Dict[str, List[Union[str, Dict[str, str]]]]]
) -> Tuple[List[pd.DataFrame], List[Union[str, Dict[str, str]]]]:
    """
    Extract and filter data from MultiQC modules based on specified section criteria.

    Args:
        module (str): The name of the MultiQC module to extract data from.
        section (Union[str, None, List, Dict]): Specification for data extraction and filtering.

    Returns:
        Tuple containing:
        - List of pandas DataFrames with extracted and filtered data
        - List of column specifications
    """
    # Retrieve all data for the specified module
    all_module_data = mqc.get_module_data(module=module)

    # Early exit if no data is found
    if not all_module_data:
        logger.warning(f"No data found for module: {module}")
        return [pd.DataFrame()], []

    # Handle simple string or None section cases
    if isinstance(section, str) or section is None:
        return extract_mqc_from_simple_section(all_module_data, section, module)

    # Handle list of strings or column specifications
    if isinstance(section, list):
        return extract_mqc_from_list_section(all_module_data, section, module)

    # Handle dictionary-based section specification
    if isinstance(section, dict):
        return extract_mqc_from_dict_section(all_module_data, section, module)

    # Fallback for unsupported section types
    logger.warning(f"Unsupported section type for module {module}: " f"type={type(section)}, value={section}")
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

    return [filter_and_rename_columns(df, columns)], columns


def extract_mqc_data(table_headers: Union[str, Path]) -> Optional[pd.DataFrame]:
    """
    Extract data from MultiQC output files.

    Args:
        mqc (object): MultiQC object.
        table_headers (Union[str, Path]): Path to the table headers file.

    Returns:
        pd.DataFrame: Extracted data
        List[Union[str, Dict[str, str]]]: List of columns containing both old and new names of columns.
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
            module_data, columns = extract_module_data(module, section)
            logger.debug("Data for %s: %s", module, module_data)

        module_data = [df.add_prefix(f"({module}) ") for df in module_data]
        columns = add_prefix_to_values_dict(columns, module)
        data.extend(module_data)
        columns_result.extend(columns)

    logger.debug("Data list: %s", data)

    return join_df(result, data) if data else result, columns_result


def write_results(contigs_mqc: pd.DataFrame, constrains_mqc: pd.DataFrame, constrains_genstats: pd.DataFrame) -> int:
    """
    Write the results to files.
    """
    samples = get_failed_samples([])
    logger.info("samples %s", samples)
    if not contigs_mqc.empty:
        logger.info("Writing Unfiltered Denovo constructs table file: contigs_overview.tsv")
        samples.extend(contigs_mqc["sample"])
        write_df(contigs_mqc.sort_values(by=["sample", "cluster", "step"]), "contigs_overview-with-iterations.tsv", [])
        table_plot = contigs_mqc[~contigs_mqc.index.isin(generate_ignore_samples(contigs_mqc))]
        write_df(table_plot.sort_values(by=["sample", "cluster", "step"]), "contigs_overview.tsv", [])

    if not constrains_mqc.empty:
        logger.info("Writing Unfiltered Mapping constructs table file: mapping_overview.tsv")
        write_df(constrains_mqc.sort_values(by=["sample", "cluster", "step"]), "mapping_overview.tsv", [])
        samples.extend(constrains_mqc["sample"])

    if not constrains_genstats.empty:
        # Add to mqc
        module = mqc.BaseMultiqcModule(name="Mapping Constrains Summary", anchor=Anchor("custom_data"))
        content = constrains_genstats.to_dict(orient="index")
        module.general_stats_addcols(content)
        mqc.report.modules.append(module)

    # Remove empty lines from the general stats data report
    samples = list(set(samples))
    mqc.report.general_stats_data = [{k: v for k, v in d.items() if k in samples} for d in mqc.report.general_stats_data]

    if mqc.report.general_stats_data:
        logger.info("Writing general stats file: samples_overview.tsv")
        samples_overview = pd.DataFrame.from_dict(get_general_stats_data_mod(), orient="index")
        samples_overview["sample"] = samples_overview.index
        write_df(reorder_columns(samples_overview, ["sample"]), "samples_overview.tsv", [])

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

    # 4. Parse our custom files into the correct tables
    overview_tables, cluster_df = load_custom_data(args)

    # 5. Make our own summary files
    # 5.1 Join with the custom contig tables
    mqc_custom_df = join_df(mqc_custom_df, overview_tables)

    if mqc_custom_df.empty:
        logger.warning("No data was found to create the contig overview table!")
        return 0

    # 5.2 reformat the dataframe
    mqc_custom_df = reformat_custom_df(mqc_custom_df, cluster_df)
    mqc_custom_df.to_csv("mqc_custom_df.after.tsv", sep="\t")

    # 5.3 split up denovo constructs and mapping (-CONSTRAIN) results
    logger.info("Splitting up denovo constructs and mapping (-CONSTRAIN) results")
    contigs_mqc, constrains_mqc = filter_constrain(mqc_custom_df, "cluster", "-CONSTRAIN")

    coalesced_constrains, constrains_genstats = reformat_constrain_df(constrains_mqc, renamed_columns, args)

    write_results(contigs_mqc, coalesced_constrains, constrains_genstats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
