#!/usr/bin/env python3

# Accept a stats file and a map file
# aggregate the stats file such that
# each category in the map file corresponds to exactly one row in the
# aggregated stats file.


import logging
import pandas as pd
import argparse
import os
import glob
import numpy as np


def usage():
    parser = argparse.ArgumentParser(description="Download from UHGG")

    parser.add_argument("-b", "--blastdir", help="directory with blast outputs")
    parser.add_argument("-o", "--output", help="final output file")
    parser.add_argument(
        "-e",
        "--ext",
        help="blast output file extension",
        default=".tsv",
        required=False,
    )
    parser.add_argument(
        "-m",
        "--mapfile",
        help="map file; tab-delimited; no header; column order ('genome_id', 'category')",
    )

    return parser.parse_args().__dict__


def parse_blast_output(filepath, extra_cols=None, extra_col_dtypes=None):
    std_blast_headers = [
        "query",
        "subj",
        "p_id",
        "aln_len",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "e_value",
        "bit_score",
        "qlen",
        "slen",
        "qcovs",
    ]
    blast_dtypes = {
        "query": str,
        "subj": str,
        "aln_len": np.float16,
        "p_id": np.float16,
        "e_value": np.float64,
        "bit_score": np.float16,
    }
    if extra_cols is not None:
        std_blast_headers.extend(extra_cols)

    if extra_col_dtypes is not None:
        blast_dtypes.update(extra_col_dtypes)

    genome_name = os.path.basename(filepath).split(".")[0]

    return pd.read_table(
        filepath,
        header=None,
        names=std_blast_headers,
        dtype=blast_dtypes,
        usecols=blast_dtypes.keys(),
    ).assign(genome_id=genome_name)


def aggregate_blast_outputs(blastdir, ext, mapfile, output):
    blast_df = pd.concat(
        (
            parse_blast_output(blastfile)
            for blastfile in glob.iglob(os.path.join(blastdir, f"*{ext}"))
        )
    ).reset_index(drop=True)

    map_df = pd.read_table(mapfile, header=None, names=["genome_id", "category"])

    summary_map_df = (
        map_df.value_counts("category")
        .reset_index()
        .rename(columns={0: "total_genomes"})
    )

    genome_blast_df = pd.merge(
        blast_df, map_df, how="left", on="genome_id"
    ).reset_index(drop=True)

    summary_blast_df = [
        summarize_results(df, min_perc_id=70).assign(category=category_name)
        for category_name, df in genome_blast_df.groupby("category")
    ]

    pd.merge(summary_blast_df, summary_map_df, on="category", how="left").to_csv(output)

    return


def summarize_results(df, min_bit_score=None, evalue_cutoff=None, min_perc_id=None):
    if min_bit_score is None:
        min_bit_score = 0

    if evalue_cutoff is None:
        evalue_cutoff = 1

    if min_perc_id is None:
        min_perc_id = 70

    return (
        df.query(
            f"(bit_score >= {min_bit_score}) and (e_value <= {evalue_cutoff}) and (p_id >= {min_perc_id})"
        )
        .groupby("query")
        .value_counts("subj")
    )


def main():
    args = usage()
    aggregate_blast_outputs(**args)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()