#!/usr/bin/env python
# Created on: Nov 21, 2023 at 1:27:12 PM
# Calculates the kmer Jaccard similarity between sequences
# About 6X faster than pairwise alignment
__author__ = "Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pyranges as pr
from mmh3 import hash
from pyfaidx import Fasta
from scipy.sparse import lil_matrix, save_npz


def seq_to_hashkmers(seq: str, k: int) -> set:
    "Converts a sequence to a set of kmers of length k"
    kmers = []
    for i in range(len(seq) - k + 1):
        kmers.append(hash(seq[i : i + k], 42, signed=False))
    return set(kmers)


def jaccard_similarity(a: set, b: set) -> float:
    "Calculates the Jaccard similarity between two sets"
    return len(a.intersection(b)) / len(a.union(b))


k = snakemake.params.k

# get transcripts from fasta
logger.info(f"Loading fasta from {snakemake.input.txome_fa}")
fasta = Fasta(snakemake.input.txome_fa, as_raw=True, sequence_always_upper=True)
txs = []
for tx in fasta.keys():
    if len(fasta[tx]) >= k:
        txs.append(tx)
    else:
        logger.warning(f"Skipping {tx}, {len(fasta[tx])} < k={k}bp")

# load gtf
logger.info(f"Loading gtf from {snakemake.input.gtf}")
t2g = (
    pr.read_gtf(snakemake.input.gtf, as_df=True)
    .query("Feature == 'transcript'")
    .set_index("transcript_id")["gene_id"]
    .to_dict()
)

# compute pairwise similarity
logger.info(f"Computing pairwise Jaccard similarity between kmers of length {k}")
matrix = lil_matrix((len(txs), len(txs)), dtype="float32")

for i, t1 in enumerate(txs):
    h1 = seq_to_hashkmers(fasta[t1][:], k)
    for j, t2 in enumerate(txs):
        if i > j:
            continue
        h2 = seq_to_hashkmers(fasta[t2][:], k)
        gene = t2g[t1] if t2g[t1] == t2g[t2] else "NA"
        score = jaccard_similarity(h1, h2)
        if score > 0:
            matrix[i, j] = score

# save matrix
logger.info(f"Saving matrix to {snakemake.output.npz}")
save_npz(snakemake.output.npz, matrix.tocsr())
