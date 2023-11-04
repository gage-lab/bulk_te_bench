#!/usr/bin/env python
# Created on: Oct 24, 2023 at 9:17:15 PM
__author__ = "Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG,
)

logger = logging.getLogger(__name__)

import shutil
from collections import defaultdict
from math import ceil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from snakemake.shell import shell


def polyester_producer(
    txome: list[SeqRecord], counts: pd.DataFrame, chunk_size: int, outdir: str
):
    """Split the counts and txome into chunks and return for polyester"""

    for i in range(0, counts.shape[0], chunk_size):
        this_counts = counts.iloc[i : i + chunk_size]

        # create a temporary fasta file with the transcripts in this chunk
        this_txome = [tx for tx in txome if tx.id in this_counts.index]
        assert (
            len(this_txome) == this_counts.shape[0]
        ), "Not all transcripts in counts were found in txome"
        SeqIO.write(this_txome, f"{outdir}/polyester_{i}.fa", "fasta")

        # convert counts to R matrix
        counts_flat = []
        for l in this_counts.transpose().to_numpy().tolist():
            for j in l:
                counts_flat.append(j)
        counts_mat = ro.r.matrix(
            counts_flat, nrow=this_counts.shape[0], ncol=this_counts.shape[1]
        )

        yield f"{outdir}/polyester_{i}.fa", counts_mat, f"{outdir}/polyester_{i}"


outdir = Path(snakemake.output[0])
shell("cp {snakemake.input.counts} {snakemake.output.counts}")
counts = pd.read_csv(snakemake.input.counts, sep="\t", index_col=0)

n_transcripts, n_samples = counts.shape
logger.info(
    f"Simulating reads from {n_transcripts} transcripts from {n_samples} samples with polyester"
)

# get full txome in memory
txome = [
    tx for tx in SeqIO.parse(snakemake.input.txome_fa, "fasta") if tx.id in counts.index
]

# run polyester in parallel
Path(snakemake.output.reads).mkdir(parents=True, exist_ok=True)
chunk_size = ceil(counts.shape[0] / snakemake.threads)
polyester = importr("polyester")
Parallel(n_jobs=snakemake.threads, backend="multiprocessing", verbose=10)(
    delayed(polyester.simulate_experiment_countmat)(
        fasta=fa,
        readmat=cmat,
        outdir=odir,
        paired=True,
        error_model="illumina5",
        gzip=True,
        seed=12,
        readlen=snakemake.params.readlen,
        strand_specific=snakemake.params.strand_specific,
    )
    for fa, cmat, odir in polyester_producer(
        txome, counts, chunk_size, snakemake.output.reads
    )
)

# merge the simulated reads
logger.info("Merging simulated reads")
samples = defaultdict(list)
for i in range(0, counts.shape[0], chunk_size):
    for f in Path(outdir / f"polyester_{i}").rglob("*.fasta.gz"):
        samples[outdir / f.name].append(f)

for outfile, infiles in samples.items():
    # concatenate all the files in the list
    with open(outfile, "wb") as outfile:
        for infile in infiles:
            with open(infile, "rb") as infile:
                outfile.write(infile.read())

# clean up
for f in Path(outdir).rglob("polyester_*.fa"):
    f.unlink()
for d in Path(outdir).rglob("polyester_*"):
    shutil.rmtree(d)
