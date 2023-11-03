#!/usr/bin/env python
# Created on: Nov 2, 2023 at 4:42:22 PM

import os
import random
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from src.txome import Txome

# change dir to this one
os.chdir(Path(__file__).parent)

OUTDIR = Path("../resources/chr22_l1hs_txome_v26")
GENOME_FA = Path("../resources/hg38.fa")
TX_GTF = Path("../resources/gencode.v26.basic.annotation.gtf.gz")
RMSK_TSV = Path("../resources/hg38.rmsk.tsv")


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
#
#   for sample in samples:
#       for gene in genes:
#        if ntx == 1:
#           DO SOMETHING
#        elif ntx == 2:
#           DO SOMETHING
#        elif ntx > 2:
#           DO SOMETHING


"""# Check out samples
metadata = pd.read_csv("../resources/GTEX/annotations_v8_GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", index_col=0, sep="\t")
samples = list(metadata.sample(frac=1, random_state=1).drop_duplicates(subset='SMTSD').index)

# create a dictionary mapping gene name to # of transcripts
gene_to_tx_df = pd.read_csv(
    "../resources/chr22_l1hs_txome/txome_t2g.tsv", sep="\t", header=None
)
# count duplicates in the gene column
gene_to_tx_values = gene_to_tx_df[1].value_counts()
# TODO: make this group by so we have dataframe with gene name as index and each transcript in list as a column -> dictionary
gene_to_tx = gene_to_tx_values.groupby()

file_path = "../resources/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
df = pd.read_csv(file_path, sep="\t", skiprows=2)
# get overlap of genes in GTEx and our Txome
genes_of_interest = set(df["name"]) & set(gene_to_tx_df[1])
df = df[df["name"].isin(genes_of_interest)]

# get random samples from each tissue
df = df.iloc[:, 0:6]  # TODO: FIX THIS TO BE ACTUAL RANDOM SAMPLES
counts = defaultdict(list)
for gene in df["name"]:
    counts["tx_id"].extend(gene_to_tx[gene])
    for sample in df:
        if gene_to_tx_df[gene] == 1:
            counts[saonmple].append(df.loc[gene, sample])
        elif gene_to_tx_df[gene] == 2:
            # do whatever
            # counts[sample].append(df.loc[gene,sample])
            pass
        elif gene_to_tx_df[gene] > 2:
            # do whatever
            # counts[sample].append(df.loc[gene,sample])
            pass
        else:
            print("ERROR: gene has no transcripts")

counts = pd.DataFrame(counts).set_index("tx_id")

txome.simulate_reads(counts, "sim_reads_4", n_jobs=32)"""
