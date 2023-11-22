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
from Bio import SeqIO
from mmh3 import hash


def seq_to_hashkmers(seq: str, k: int) -> set:
    "Converts a sequence to a set of kmers of length k"
    kmers = []
    for i in range(len(seq) - k + 1):
        kmers.append(hash(seq[i : i + k], 42, signed=False))
        # kmers.append(seq[i:i+k])
    return set(kmers)


def jaccard_similarity(a: set, b: set) -> float:
    "Calculates the Jaccard similarity between two sets"
    return len(a.intersection(b)) / len(a.union(b))


# How big does this get with the full txome?
# assuming all txs are 6kb (conservative, most are shorter)
# each tx will have 6000-k+1 kmers, hashed to 32 bits (4 bytes) each
# 6000 * 4 = 24kb per tx
# 24kb * (roughly) 100,000 txs = 2.4GB
logger.info("Hashing kmers from transcripts")
hashed = {}
for r in SeqIO.parse(snakemake.input.txome_fa, "fasta"):
    hashed[r.id] = seq_to_hashkmers(str(r.seq).upper(), k=50)

# load gtf
t2g = (
    pr.read_gtf(snakemake.input.gtf, as_df=True)
    .query("Feature == 'transcript'")
    .set_index("transcript_id")["gene_id"]
    .to_dict()
)

logger.info("Computing pairwise Jaccard similarity")
with open(snakemake.output.tsv, "w") as f:
    f.write("tx1\ttx2\tgene\tjaccard_similarity\n")
    for i, t1 in enumerate(hashed):
        for j, t2 in enumerate(hashed):
            if i > j:
                continue
            gene = t2g[t1] if t2g[t1] == t2g[t2] else "NA"
            score = jaccard_similarity(hashed[t1], hashed[t2])
            f.write(f"{t1}\t{t2}\t{gene}\t{score}\n")
