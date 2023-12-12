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
import pandas as pd
import pyranges as pr
from Bio import SeqIO
from mmh3 import hash
from scipy.stats import dirichlet


# TODO: define function to make tpm from counts
# TODO: simulate splicing rate per gene
def tx_from_gene_counts(gene_counts: pd.DataFrame, g2t: dict, tx: list[str]):
    """
    Generate transcript counts from gene counts using dirichlet distribution to split counts between isoforms.
    :param gene_counts: pd.DataFrame with gene counts
    :param g2t: dict mapping gene to transcripts
    :param tx: list of all transcripts
    """

    counts = pd.DataFrame(0, index=tx, columns=gene_counts.columns)
    for gene in gene_counts.index:
        ntx = len(g2t[gene])
        for sample in gene_counts.columns:
            # set unique seed for each gene/sample pair
            seed = hash(gene + sample, 42, signed=False)
            np.random.seed(seed)
            gc = gene_counts.loc[gene, sample]
            cointoss = np.random.randint(0, 2)  # 0 or 1
            if ntx == 1:
                tx_counts = np.array([gc])
            elif ntx == 2:
                # For genes with 2 transcripts, gene counts are either
                # (i) split among both isoforms according to flat Dirichlet distribution (α = (1,1))
                # (ii) attributed to a single isoform.
                if cointoss == 0:
                    tx_counts = dirichlet.rvs([1, 1], random_state=seed)[0] * gc
                else:
                    tx_counts = np.array([gc, 0])

            elif ntx > 2:
                # For genes with more than 2 transcripts, gene counts are either
                # (i) split among 3 randomly chosen isoforms according to flat Dirichlet distribution (α = (1,1,1))
                # (ii) attributed to a single isoform.
                if cointoss == 0:
                    tx_counts = dirichlet.rvs([1, 1, 1], random_state=seed)[0] * gc
                else:
                    tx_counts = np.array([gc, 0, 0])

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

    return counts


# read in gtfs
genes_gtf = pr.read_gtf(snakemake.input.genes_gtf, as_df=True).query(
    "Feature == 'transcript'"
)
# genes_gtf = genes_gtf[~genes_gtf["transcript_id"].str.contains("-I")] # exclude unspliced transcripts
g2t = (
    genes_gtf.groupby("gene_id").apply(lambda x: x["transcript_id"].tolist()).to_dict()
)
rmsk = pr.read_gtf(snakemake.input.rmsk_gtf, as_df=True)

# setup variables
txs = genes_gtf["transcript_id"].unique()
samples = [f"sample_{i:02d}" for i in range(1, int(snakemake.params.nsamples) + 1)]

## UNIFORM COUNTS ##
if snakemake.wildcards.txte_sim == "uniform_sim":

    logger.info("Generating uniform counts")
    counts = pd.DataFrame(0, index=txs, columns=samples)
    for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
        if tx.id not in txs:
            continue
        for sample in samples:
            counts.loc[tx.id, sample] = 3 * len(tx.seq) // 100

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

    logger.info(
        f"{len(shared_genes)}/{len(gtex_counts)} GTEx genes are shared with txome\n"
        f"{len(shared_genes)}/{len(g2t)} txome genes are shared with GTEx"
    )

    # subset GTEx to shared genes and samples
    gtex_counts = gtex_counts.set_index("Name").loc[list(shared_genes), samples]
    logger.info(f"GTEx counts shape after subsetting: {gtex_counts.shape}")

    # generate the transcript level counts from GTEx gene counts
    logger.info("Generating transcript level counts from GTEx gene counts")
    counts = tx_from_gene_counts(gtex_counts, g2t, all_tx)

    # only keep desired number of samples
    counts = counts.iloc[:, 0 : int(snakemake.params.nsamples)]

elif snakemake.wildcards.txte_sim == "single_intergenic_l1hs":

    # filter for l1hs
    rmsk = rmsk[rmsk["gene_id"] == "L1HS"]
    tes = rmsk["transcript_id"].unique()

    # error check
    if rmsk.empty:
        logger.error("The txome does not have any L1HS loci.")
        raise ValueError("The txome does not have any L1HS loci.")
    if rmsk["contained_in"].isna().sum() < 1:
        logger.error("The txome does not have any Intergenic TE loci.")
        raise ValueError("The txome does not have any Intergenic TE loci.")

    # for each sample, randomly choose a l1hs intergenic locus to add counts to
    logger.info("Simulating single intergenic l1hs locus expression for each sample")
    counts = pd.DataFrame(0, index=tes, columns=samples)
    for sample in samples:
        my_te = np.random.choice(rmsk[rmsk["contained_in"].isna()].transcript_id)
        for te in rmsk.itertuples():
            if te.transcript_id == my_te:
                # avg tpm in GTEx per tx in 40-120
                counts.loc[te.transcript_id, sample] = (
                    20 * np.abs(te.End - te.Start) // 100
                )

else:
    logger.error(
        f"The simulation name {snakemake.wildcards.txte_sims} is not supported."
    )
    raise ValueError(
        f"The simulation name {snakemake.wildcards.txte_sims} is not supported."
    )

# rename columns to sample_01, sample_02, etc.
counts.to_csv(snakemake.output.counts, sep="\t", index=True)
