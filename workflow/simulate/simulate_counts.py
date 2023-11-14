#!/usr/bin/env python
# Created on: Nov 8, 2023 at 3:01:15 PM
__author__ = ["Michael Cuoco", "Joelle Fabyishenko"]

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

# read in gtfs
genes_gtf = pr.read_gtf(snakemake.input.genes_gtf)
rmsk_gtf = pr.read_gtf(snakemake.input.rmsk_gtf)

# label rmsk elements as intergenic or intronic
rmsk = rmsk_gtf.join(genes_gtf, how="left", report_overlap=True, suffix="_gene").df
rmsk["Intronic"] = rmsk["Overlap"] == (rmsk["End"] - rmsk["Start"])
rmsk["Intergenic"] = rmsk["Overlap"] < 0

# convert to df
rmsk_gtf = rmsk_gtf.df
genes_gtf = genes_gtf.df

## TEST COUNTS ##
if snakemake.wildcards.txte_sim == "test_sim":
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
        for sample in range(0, int(snakemake.params.nsamples)):
            counts[sample].append(20 * len(tx.seq) // 100)

    counts = pd.DataFrame(counts)

## UNIFORM COUNTS ##
elif snakemake.wildcards.txte_sim == "uniform_sim":

    logger.info("Generating uniform counts")
    counts = defaultdict(list)
    for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
        if tx.id not in genes_gtf["transcript_id"].unique():
            continue
        counts["tx_id"].append(tx.id)
        for sample in range(0, int(snakemake.params.nsamples)):
            counts[sample].append(20 * len(tx.seq) // 100)

    counts = pd.DataFrame(counts)

## GTEx COUNTS ##
elif snakemake.wildcards.txte_sim == "gtex_sim":
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
    counts = pd.DataFrame(0, index=all_tx, columns=gtex_counts.columns)
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
    # only keep desired number of samples
    counts = counts.iloc[:, 0 : int(snakemake.params.nsamples)].reset_index()


elif snakemake.wildcards.txte_sim == "single_intergenic_l1hs":

    # filter for l1hs
    rmsk = rmsk[rmsk["gene_id"] == "L1HS"]

    # error check
    if rmsk.empty:
        logger.error("The txome does not have any L1HS loci.")
        raise ValueError("The txome does not have any L1HS loci.")
    if rmsk["Intergenic"].sum() < 1:
        logger.error("The txome does not have any Intergenic TE loci.")
        raise ValueError("The txome does not have any Intergenic TE loci.")

    # for each sample, randomly choose a l1hs intergenic locus to add counts to
    logger.info("Simulating single intergenic l1hs locus expression for each sample")
    counts = pd.DataFrame(
        0,
        index=rmsk.transcript_id.unique(),
        columns=[i for i in range(0, int(snakemake.params.nsamples))],
    )
    for sample in counts.columns:
        my_te = np.random.choice(rmsk[rmsk["Intergenic"]].transcript_id)
        for te in rmsk.itertuples():
            if te.transcript_id == my_te:
                counts.loc[te.transcript_id, sample] = (
                    20 * np.abs(te.End - te.Start) // 100
                )
    counts = counts.reset_index().rename(columns={"transcript_id": "tx_id"})

else:
    logger.error(
        f"The simulation name {snakemake.wildcards.txte_sims} is not supported."
    )
    raise ValueError(
        f"The simulation name {snakemake.wildcards.txte_sims} is not supported."
    )

counts.to_csv(snakemake.output.counts, sep="\t", index=False)
