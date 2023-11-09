#!/usr/bin/env python
# Created on: Nov 8, 2023 at 5:19:15 PM
__author__ = "Joelle Faybishenko, Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import numpy as np

np.random.seed(1)
from collections import defaultdict

import pandas as pd
import pyranges as pr
from Bio import SeqIO

# Get all l1hs loci from rmsk
genes_gtf = pr.read_gtf(snakemake.input.genes_gtf)
rmsk_gtf = pr.read_gtf(snakemake.input.rmsk_gtf)

# Find ones that are intergenic, by seeing if there are overlaps between location of loci in rmsk and gene gtf
df = rmsk_gtf.join(genes_gtf, how="left", report_overlap=True, suffix="_rmsk").df
df["Intergenic"] = df["Overlap"] == (df["End"] - df["Start"])
df["Intronic"] = df["Overlap"] < 0

if df["Intergenic"].sum() < 1:
    logger.error("The txome does not have any Intergenic loci.")
    raise ValueError("The txome does not have any Intergenic loci.")
if df["Intronic"].sum() < 1:
    logger.error("The txome does not have any intronic loci.")
    raise ValueError("The txome does not have any intronic loci.")

counts_mtx = pd.read_csv(snakemake.input.counts, sep="\t")


row_values = []
tes = []
logger.info("Iterating over fasta txs")
# Get average tpm in the tx count table (inverse equation to uniform/test sims)
for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
    if tx.id in list(counts_mtx["tx_id"]):
        # Inverse equation to uniform/test sims
        row_values.append(
            np.mean(
                counts_mtx.loc[counts_mtx["tx_id"] == tx.id,].values[
                    0
                ][1:]
            )
            * 100
            // len(tx.seq)
        )
    else:
        if snakemake.wildcards.tx_sim == "test_sim":
            if "ENS" in tx.id:
                continue

        tes.append(tx)


## SIMULATION 1: single intergenic loci ##
if snakemake.wildcards.te_sim == "single_intergenic_l1hs":
    logger.info("Iterating over L1HS loci")

    te_counts = defaultdict(list)
    keep_empty = False
    # warning: this will not work for test data
    for tx in tes:
        te_counts["tx_id"].append(tx.id)
        # if intergenic locus added
        if keep_empty:
            # Add other l1hs loci (0 counts)
            for sample in counts_mtx.columns[1:]:
                te_counts[sample].append(0)

        # if first intergenic locus
        elif tx.id in list(
            df.loc[
                df["Intergenic"],
            ]["transcript_id"]
        ):
            # Add row to table with 1 intergenic l1hs that has the normalized average counts for all samples
            for sample in counts_mtx.columns[1:]:
                te_counts[sample].append(np.mean(row_values) * len(tx.seq) // 100)
            keep_empty = True

        # if not first intergenic locus nor intergenic locus added
        else:
            # Add other l1hs loci (0 counts)
            for sample in counts_mtx.columns[1:]:
                te_counts[sample].append(0)

else:
    pass

pd.concat([counts_mtx, pd.DataFrame(te_counts)], axis=0).to_csv(
    snakemake.output.te_counts, sep="\t", index=False
)
