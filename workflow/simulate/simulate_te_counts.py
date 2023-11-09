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

# Get mean tpm in the tx count table
counts = pd.read_csv(snakemake.input.tx_counts, sep="\t", index_col=0)
mean_counts = counts.mean(axis=1)

tx_tpms = []
logger.info("Iterating over fasta txs")
for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
    if tx.id in list(counts.index):
        # Calculate tpm from counts
        tx_tpms.append(mean_counts[tx.id] * 100 // len(tx.seq))

mean_tx_tpm = np.mean(tx_tpms)

# Get all l1hs loci from rmsk
genes_gtf = pr.read_gtf(snakemake.input.genes_gtf)
rmsk_gtf = pr.read_gtf(snakemake.input.rmsk_gtf)

# Find ones that are intergenic, by seeing if there are overlaps between location of loci in rmsk and gene gtf
rmsk = rmsk_gtf.join(genes_gtf, how="left", report_overlap=True, suffix="_gene").df
rmsk["Intronic"] = rmsk["Overlap"] == (rmsk["End"] - rmsk["Start"])
rmsk["Intergenic"] = rmsk["Overlap"] < 0


def contains_intergenic_intronic(df):
    if df["Intergenic"].sum() < 1:
        logger.error("The txome does not have any Intergenic TE loci.")
        raise ValueError("The txome does not have any Intergenic TE loci.")
    if df["Intronic"].sum() < 1:
        logger.error("The txome does not have any intronic TE loci.")
        raise ValueError("The txome does not have any intronic TE loci.")


contains_intergenic_intronic(rmsk)

## SIMULATION 1: single intergenic loci ##
if snakemake.wildcards.te_sim == "single_intergenic_l1hs":
    logger.info("Iterating over L1HS loci")

    # filter for l1hs
    rmsk = rmsk[rmsk["gene_id"] == "L1HS"]
    if rmsk.empty:
        logger.error("The txome does not have any L1HS loci.")
        raise ValueError("The txome does not have any L1HS loci.")
    contains_intergenic_intronic(rmsk)

    # for each sample, randomly choose a l1hs intergenic locus to add counts to
    for sample in counts.columns:
        my_te = np.random.choice(rmsk[rmsk["Intergenic"]].transcript_id)
        for te in rmsk.itertuples():
            if te.transcript_id == my_te:
                counts.loc[te.transcript_id, sample] = (
                    mean_tx_tpm * np.abs(te.End - te.Start) // 100
                )
            else:
                counts.loc[te.transcript_id, sample] = 0
else:
    logger.error("The simulation type is not supported.")
    raise ValueError("The simulation type is not supported.")

counts.to_csv(snakemake.output.counts, sep="\t")
