#!/usr/bin/env python
# Created on: Nov 13, 2023 at 3:03:51 PM
__author__ = "Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import gzip
import random
from pathlib import Path

from Bio import SeqIO

random.seed(1)
from itertools import chain

import pandas as pd
import pyranges as pr

# concat count matrices
logger.info("Concatenating count matrices")
t2g = (
    pr.read_gtf(snakemake.input.gtf, as_df=True)
    .query("Feature == 'transcript'")
    .rename(columns={"transcript_id": "tx_id"})
    .set_index("tx_id")["gene_id"]
    .to_dict()
)
te_counts = pd.read_csv(snakemake.input.te_counts[0], sep="\t", index_col=0)
tx_counts = pd.read_csv(snakemake.input.tx_counts[0], sep="\t", index_col=0)
counts = pd.concat([te_counts, tx_counts], axis=0)
counts.to_csv(snakemake.output.tx_counts, sep="\t")

# add a gene id column
counts["gene_id"] = counts.index.map(t2g)
counts.groupby("gene_id").sum().to_csv(snakemake.output.gene_counts, sep="\t")

# concat reads
logger.info("Concatenating tx and te reads")

te_r1_files = sorted(Path(snakemake.input.te_reads[0]).rglob("*_1.fasta.gz"))
te_r2_files = sorted(Path(snakemake.input.te_reads[0]).rglob("*_2.fasta.gz"))
tx_r1_files = sorted(Path(snakemake.input.tx_reads[0]).rglob("*_1.fasta.gz"))
tx_r2_files = sorted(Path(snakemake.input.tx_reads[0]).rglob("*_2.fasta.gz"))

Path(snakemake.output.reads).mkdir(parents=True, exist_ok=True)
for files in [(te_r1_files, tx_r1_files), (te_r2_files, tx_r2_files)]:
    assert len(files[0]) == len(files[1]), "Number of read files must be equal"
    for (
        tx,
        te,
    ) in zip(*files):
        assert tx.name == te.name, "Read files must be in the same order"
        outfile = Path(snakemake.output.reads) / tx.name
        outfile.touch()
        logger.info(f"Concatenating {tx.name}")

        # read and shuffle the reads
        tx_in = gzip.open(tx, "rt")
        tx_reads = SeqIO.parse(tx_in, "fasta")
        te_in = gzip.open(te, "rt")
        te_reads = SeqIO.parse(te_in, "fasta")
        reads = chain(tx_reads, te_reads)

        # write the reads
        with gzip.open(outfile, "wt") as r_out:
            for read in sorted(reads, key=lambda x: random.random()):
                SeqIO.write(read, r_out, "fasta")

        # close the files
        tx_in.close()
        te_in.close()
