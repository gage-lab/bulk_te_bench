#!/usr/bin/env python
# Created on: Nov 2, 2023 at 4:42:22 PM

import os
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import dirichlet

from src.txome import Txome

random.seed(1)
np.random.seed(1)


def make_realistic_counts(
    GTEx_counts: Path(),
    GTEx_metadata: Path(),
    GENCODE_annotation: Path(),
    outdir=Path(),
    txome_exists: bool = False,
    chromosone: str = "chr22",
) -> pd.DataFrame():
    """
    Simulate realistic counts from GTEx data
    @param GTEx_counts: Path to GTEx counts file
    @param GTEx_metadata: Path to GTEx metadata file
    @param GENCODE_annotation: Path to GENCODE annotation file
    @param outdir: Path to transcriptome output directory
    @param txome_exists: Boolean indicating whether or not the transcriptome already exists
    @param chromosone: Chromosome to subset transcriptome
    @return simulated transcript count matrix
    """

    # change dir to this one
    os.chdir(Path(__file__).parent)

    OUTDIR = outdir
    GENOME_FA = Path("../resources/hg38.fa")
    RMSK_TSV = Path("../resources/hg38.rmsk.tsv")
    TX_GTF = GENCODE_annotation
    GTEx_COUNTS_PATH = GTEx_counts
    GTEx_METADATA_PATH = GTEx_metadata

    txome = Txome(
        outdir=OUTDIR, genome_fa=GENOME_FA, tx_gtf=TX_GTF, chromosome=chromosone
    )
    txome.get_rmsk_seqs(
        parsed_rmsk=RMSK_TSV,
        query="repName == 'L1HS' & length > 6000",
    )

    txome.txome_fa = OUTDIR / "txome.fa"

    if not txome_exists:
        txome.make_txome()

    all_tx = set([tx.id for tx in SeqIO.parse(txome.txome_fa, "fasta")])

    # create a dictionary mapping gene name to # of transcripts
    gene_to_tx_df = pd.read_csv(OUTDIR / "txome_t2g.tsv", sep="\t", header=None)
    # subset to only transcripts in txome
    gene_to_tx_df = gene_to_tx_df.loc[gene_to_tx_df[0].isin(all_tx), :]
    # Mapping of genes to list of transcripts
    gene_to_tx = gene_to_tx_df.groupby(1).apply(lambda x: x[0].tolist()).to_dict()

    full_GTEx = pd.read_csv(GTEx_COUNTS_PATH, sep="\t", skiprows=2)

    # Randomly sample from each tissue
    metadata = pd.read_csv(
        GTEx_METADATA_PATH,
        index_col=0,
        sep="\t",
    )
    metadata = metadata.loc[(list(full_GTEx.columns)[2:]), "SMTS"]
    samples = list(metadata.sample(frac=1, random_state=1).drop_duplicates().index)

    # Get overlaps of genes in GTEx and our Txome
    genes_of_interest = set(full_GTEx["Name"]) & set(gene_to_tx_df[1])
    full_GTEx.set_index("Name", inplace=True)
    subset_GTEx = full_GTEx.loc[full_GTEx.index.isin(genes_of_interest), samples]

    # Voodoo: get realistic tx counts from GTEx genes
    counts = defaultdict(list)
    for gene in subset_GTEx.index:
        counts["tx_id"].extend(gene_to_tx[gene])
        for sample in subset_GTEx.columns:
            if len(gene_to_tx[gene]) == 1:
                counts[sample].append(subset_GTEx.loc[gene, sample])
            elif len(gene_to_tx[gene]) == 2:
                TPM = subset_GTEx.loc[gene, sample]

                cointoss = random.randint(0, 1)
                if cointoss == 0:
                    # use dirichlet to split counts between two isoforms
                    distribution = (
                        dirichlet.rvs([1, 1], size=1)[0] * TPM
                    )  # multiply this by counts
                else:
                    # give to one isoform
                    distribution = [TPM, 0]

                np.random.shuffle(distribution)
                counts[sample].extend(distribution)

            elif len(gene_to_tx[gene]) > 2:
                # TPMs were either (i) split among three randomly chosen isoforms according to a flat Dirichlet distribution
                # (α = (1,1,1)) or (ii) attributed to a single isoform.
                TPM = subset_GTEx.loc[gene, sample]

                cointoss = random.randint(0, 1)
                if cointoss == 0:
                    # use dirichlet to split counts between two isoforms
                    distribution = np.array(
                        dirichlet.rvs([1, 1, 1], size=1)[0] * TPM
                    )  # multiply this by counts
                else:
                    # give to one isoform
                    distribution = np.array([TPM, 0, 0])

                # if there is more than 3 transcripts, add zeros to the distribution
                if len(distribution) != len(gene_to_tx[gene]):
                    distribution = np.pad(
                        distribution, (0, (len(gene_to_tx[gene]) - len(distribution)))
                    )

                # randomize
                np.random.shuffle(distribution)
                counts[sample].extend(distribution)
            else:
                print("ERROR: gene has no transcripts")

    # get non-duplicated transcripts in txome
    all_tx = set([tx.id for tx in SeqIO.parse(txome.txome_fa, "fasta")])
    counts = pd.DataFrame(counts).set_index("tx_id")

    return counts
