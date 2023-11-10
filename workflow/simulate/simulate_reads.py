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

import gzip
import shutil
from collections import defaultdict
from math import ceil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()


def polyester_producer(
    txome: list[SeqRecord], counts: pd.DataFrame, chunk_size: int, outdir: str
):
    """Split the counts and txome into chunks and return for polyester"""

    for i in range(0, counts.shape[0], chunk_size):
        this_counts = counts.iloc[i : i + chunk_size]

        # create a temporary fasta file with the transcripts in this chunk
        this_txome = []
        for ctx in this_counts.index:
            for ttx in txome:
                if ttx.id == ctx:
                    this_txome.append(ttx)
        SeqIO.write(this_txome, f"{outdir}/polyester_{i}.fa", "fasta")

        # convert counts to R matrix
        counts_mat = pandas2ri.py2rpy(this_counts)
        counts_mat = ro.r["as.matrix"](counts_mat)

        yield f"{outdir}/polyester_{i}.fa", counts_mat, f"{outdir}/polyester_{i}"


outdir = Path(snakemake.output[0])
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


# check that the number of reads in the output is correct
logger.info("Checking that the number of reads in the output is correct")
read_counts = pd.DataFrame(0, index=counts.index, columns=counts.columns, dtype=float)
for fagz in samples.keys():

    # ignore read 2
    if "_2.fasta.gz" in str(fagz):
        continue

    # get sample name
    sample = str(int(fagz.name[7:9]) - 1)

    # count reads from each tx
    with gzip.open(fagz, "rt") as fa:
        for read in SeqIO.parse(fa, "fasta"):
            tx = read.id.split("/")[1].split(";")[0]
            read_counts.loc[tx, sample] += 1

if not read_counts.equals(counts):
    logger.error("Simulated reads do not match simulate counts")
    raise ValueError("Simulated reads do not match simulate counts")

# clean up
for f in Path(outdir).rglob("polyester_*.fa"):
    f.unlink()
for d in Path(outdir).rglob("polyester_*"):
    shutil.rmtree(d)
