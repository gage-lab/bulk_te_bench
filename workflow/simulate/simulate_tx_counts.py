#!/usr/bin/env python
# Created on: Nov 8, 2023 at 3:01:15 PM
__author__ = "Michael Cuoco"

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
from scipy.stats import dirichlet

genes_gtf = pr.read_gtf(snakemake.input.genes_gtf, as_df=True)
rmsk_gtf = pr.read_gtf(snakemake.input.rmsk_gtf, as_df=True)

## TEST COUNTS ##
if snakemake.wildcards.tx_sim == "test_sim":
    logger.info("Generating test counts")
    TX = [
        "ENST00000401850.5",
        "ENST00000401959.6",
        "ENST00000401975.5",
        "ENST00000401994.5",
    ]

    counts = defaultdict(list)
    for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
        if tx.id not in TX:
            continue
        counts["tx_id"].append(tx.id)
        for sample in range(0, 2):
            counts[sample].append(20 * len(tx.seq) // 100)

    pd.DataFrame(counts).to_csv(snakemake.output.tx_counts, sep="\t", index=False)

## UNIFORM COUNTS ##
elif snakemake.wildcards.tx_sim == "uniform_sim":
    logger.info("Generating uniform counts")
    counts = defaultdict(list)
    for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
        if tx.id not in genes_gtf["transcript_id"].unique():
            continue
        counts["tx_id"].append(tx.id)
        for sample in range(0, 7):
            counts[sample].append(20 * len(tx.seq) // 100)
    pd.DataFrame(counts).to_csv(snakemake.output.tx_counts, sep="\t", index=False)


## GTEx COUNTS ##
elif snakemake.wildcards.tx_sim == "gtex_sim":
    logger.info("Generating GTEx counts")

    # take one sample from each tissue
    # takes ~5 min to read in the gtex counts
    gtex_counts = pd.read_csv(snakemake.input.gtex_counts, sep="\t", skiprows=2)
    logger.info(f"GTEx counts shape: {gtex_counts.shape}")
    samples = (
        pd.read_csv(snakemake.input.gtex_metadata, sep="\t", index_col=0)
        .loc[gtex_counts.columns[2:], "SMTS"]
        .sample(frac=1, random_state=1)
        .drop_duplicates()
        .index
    )
    logger.info(f"Grabbing 1 GTEx sample from each tissue: {samples} samples")

    # get shared genes in GTEx and our txome
    shared_genes = set(gtex_counts["Name"]) & set(genes_gtf["gene_id"])
    all_tx = genes_gtf.query("Feature == 'transcript'")["transcript_id"].unique()
    genes_gtf_only_genes = set(genes_gtf["gene_id"]) - set(gtex_counts["Name"])

    g2t = (
        genes_gtf.query("Feature == 'transcript'")
        .groupby("gene_id")
        .apply(lambda x: x["transcript_id"].tolist())
        .to_dict()
    )

    logger.info(
        f"{len(shared_genes)}/{len(gtex_counts)} GTEx genes are shared with txome"
    )
    logger.info(f"{len(shared_genes)}/{len(g2t)} txome genes are shared with GTEx")

    # subset GTEx to shared genes and samples
    gtex_counts = gtex_counts.set_index("Name").loc[list(shared_genes), samples]
    logger.info(f"GTEx counts shape after subsetting: {gtex_counts.shape}")

    # generate the transcript level counts from GTEx gene counts
    logger.info("Generating transcript level counts from GTEx gene counts")
    counts = pd.DataFrame(index=all_tx, columns=gtex_counts.columns)
    for gene in gtex_counts.index:
        ntx = len(g2t[gene])
        for sample in gtex_counts.columns:
            gene_counts = gtex_counts.loc[gene, sample]
            if ntx == 1:
                tx_counts = np.array([gene_counts])
            elif ntx == 2:
                cointoss = np.random.randint(0, 2)  # 0 or 1
                if cointoss == 0:
                    # use dirichlet to split counts between two isoforms
                    tx_counts = (
                        dirichlet.rvs([1, 1], size=1)[0] * gene_counts
                    )  # multiply this by counts
                else:
                    # give to one isoform
                    tx_counts = np.array([gene_counts, 0])

            elif ntx > 2:
                # Gene counts either
                # (i) split among 3 randomly chosen isoforms according to flat Dirichlet distribution (α = (1,1,1))
                # (ii) attributed to a single isoform.

                cointoss = np.random.randint(0, 2)
                if cointoss == 0:
                    # use dirichlet to split counts two isoforms
                    tx_counts = np.array(
                        dirichlet.rvs([1, 1, 1], size=1)[0] * gene_counts
                    )  # multiply this by counts
                else:
                    # give to one isoform
                    tx_counts = np.array([gene_counts, 0, 0])

                # if there is more than 3 transcripts, add zeros to the tx_counts
                if len(tx_counts) != ntx:
                    tx_counts = np.pad(tx_counts, (0, (ntx - len(tx_counts))))
                    if len(tx_counts) != ntx:
                        logger.error("ERROR: tx_counts != ntx")
                        raise ValueError("ERROR: tx_counts != ntx")

            else:
                raise ValueError("ERROR: gene has no transcripts")

            np.random.shuffle(tx_counts)
            counts.loc[g2t[gene], sample] = tx_counts.astype(int)

    # rename columns to numbers
    counts.columns = [i for i in range(0, counts.shape[1])]

    # for any genes that are in our txome but not in GTEx, set counts to 0
    counts.fillna(0).to_csv(snakemake.output.tx_counts, sep="\t")

else:
    logger.error("Unknown simulation name!")
    raise ValueError("Unknown simulation name!")
