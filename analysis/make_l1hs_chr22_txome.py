#!/usr/bin/env python
# Created on: Oct 24, 2023 at 9:17:15 PM
__author__ = "Michael Cuoco"

import os
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from src.txome import Txome

# change dir to this one
os.chdir(Path(__file__).parent)

OUTDIR = Path("../resources/chr22_l1hs_txome")
GENOME_FA = Path("../resources/hg38.fa")
TX_GTF = Path("../resources/gencode.v44.primary_assembly.basic.annotation.gtf")
RMSK_TSV = Path("../resources/hg38.rmsk.tsv")


txome = Txome(outdir=OUTDIR, genome_fa=GENOME_FA, tx_gtf=TX_GTF, chromosome="chr22")
txome.get_rmsk_seqs(
    parsed_rmsk=RMSK_TSV,
    query="repName == 'L1HS' & length > 6000",
)
txome.txome_fa = Path("../resources/chr22_l1hs_txome/txome.fa")
txome.make_txome()

# SIMULATION 1: uniform FPKM, hold all transcripts constant across samples, only increase L1HS expression
# make count matrix to simulate from
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
counts = defaultdict(list)
for tx in SeqIO.parse(txome.txome_fa, "fasta"):
    counts["tx_id"].append(tx.id)
    for sample in range(0, 6):
        if "ENS" in tx.id:
            counts[sample].append(20 * len(tx.seq) // 100)
        elif "chr" in tx.id:
            counts[sample].append(sample * len(tx.seq) // 100)
        else:
            counts[sample].append(0)

counts = pd.DataFrame(counts).set_index("tx_id")

txome.simulate_reads(counts, "sim_reads_1", n_jobs=32)

# TODO: add additional simulations here
