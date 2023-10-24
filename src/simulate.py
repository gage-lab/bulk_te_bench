#!/usr/bin/env python
# Created on: Oct 24, 2023 at 1:32:22 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed
from math import ceil
from collections import defaultdict
import shutil
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr


def producer(txome: list, counts: pd.DataFrame, chunk_size: int, outdir: str):
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


def run_polyester(txome_fa: str, counts: pd.DataFrame, n_jobs: int, outdir: str):
    """
    Run polyester to simulate reads from a transcriptome in parallel
    :param txome_fa: path to transcriptome fasta
    :param counts: pandas dataframe with reads per transcript per sample
    :param n_jobs: number of parallel threads to use
    :param outdir: output directory
    """
    # check input
    if not Path(txome_fa).exists():
        raise FileNotFoundError(f"No such file {txome_fa}")

    # make output directory
    Path(outdir).mkdir(exist_ok=True)
    counts.to_csv(f"{outdir}/counts.tsv", sep="\t")

    n_transcripts, n_samples = counts.shape
    logger.info(
        f"Simulating reads from {n_transcripts} transcripts from {n_samples} samples with polyester"
    )

    # get full txome in memory
    txome = [tx for tx in SeqIO.parse(txome_fa, "fasta") if tx.id in counts.index]

    # run polyester in parallel
    chunk_size = ceil(counts.shape[0] / n_jobs)
    polyester = importr("polyester")
    Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=10)(
        delayed(polyester.simulate_experiment_countmat)(
            fasta=fa,
            readmat=cmat,
            outdir=odir,
            paired=True,
            error_model="illumina5",
            gzip=True,
            seed=12,
        )
        for fa, cmat, odir in producer(txome, counts, chunk_size, outdir)
    )

    # merge the simulated reads
    logger.info("Merging simulated reads")
    samples = defaultdict(list)
    for i in range(0, counts.shape[0], chunk_size):
        for f in Path(f"{outdir}/polyester_{i}").rglob("*.fasta.gz"):
            samples[f"{outdir}/{f.name}"].append(f)

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
