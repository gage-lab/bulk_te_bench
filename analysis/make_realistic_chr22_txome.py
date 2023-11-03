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

# change dir to this one
os.chdir(Path(__file__).parent)

random.seed(1)
np.random.seed(1)

OUTDIR = Path("../resources/chr22_l1hs_txome_v26")
GENOME_FA = Path("../resources/hg38.fa")
TX_GTF = Path("../resources/gencode.v26.adjusted.basic.annotation.chr22.gtf.gz")
RMSK_TSV = Path("../resources/hg38.rmsk.tsv")
GTEx_COUNTS_PATH = Path(
    "../resources/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
)
GTEx_METADATA_PATH = Path(
    "../resources/GTEX/annotations_v8_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
)


txome = Txome(outdir=OUTDIR, genome_fa=GENOME_FA, tx_gtf=TX_GTF, chromosome="chr22")
txome.get_rmsk_seqs(
    parsed_rmsk=RMSK_TSV,
    query="repName == 'L1HS' & length > 6000",
)
txome.txome_fa = Path("../resources/chr22_l1hs_txome_v26/txome.fa")
txome.make_txome()


#   SIMULATION 5: realistic counts
#   Simulate reads with realistic counts
#   1. read in GTEX gene count matrix and metadata
#   2. Select transcripts on chr22 only
#   3. Intersect gene IDs from GTEx quantification to our Txome's t2g.tsv file (or switch to their version of GENCODE and remake Txome)
#   4. Take 1 sample from each tissue (ensure it's the same samples each time this script is run)
#   5. Add L1HS transcripts to the count matrix with a constant read count (maybe median counts for that sample?)


# create a dictionary mapping gene name to # of transcripts
gene_to_tx_df = pd.read_csv(
    "../resources/chr22_l1hs_txome_v26/txome_t2g.tsv", sep="\t", header=None
)
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
chr22_GTEx = full_GTEx.loc[full_GTEx.index.isin(genes_of_interest), samples]


# Voodoo: get realistic tx counts from GTEx genes
counts = defaultdict(list)
for gene in chr22_GTEx.index:
    counts["tx_id"].extend(gene_to_tx[gene])
    for sample in chr22_GTEx.columns:
        if len(gene_to_tx[gene]) == 1:
            counts[sample].append(chr22_GTEx.loc[gene, sample])
        elif len(gene_to_tx[gene]) == 2:
            transcripts = gene_to_tx[gene]
            TPM = chr22_GTEx.loc[gene, sample]

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
            transcripts = gene_to_tx[gene]
            TPM = chr22_GTEx.loc[gene, sample]

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


# save and simulate
counts = pd.DataFrame(counts).set_index("tx_id")

txome.simulate_reads(counts, "sim_reads_5", n_jobs=32)
