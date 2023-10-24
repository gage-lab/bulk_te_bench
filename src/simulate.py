from pathlib import Path
import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed
from math import ceil
from collections import defaultdict
import logging, shutil


def _polyester(txome_fa: str, counts: pd.DataFrame, outdir: str) -> None:
    """
    Run polyester to simulate reads from a transcriptome
    :param txome_fa: path to transcriptome fasta
    :param counts: pandas dataframe with transcript counts
    :param outdir: output directory
    :return: None
    """
    from rpy2 import robjects as ro
    from rpy2.robjects.packages import importr

    # make output directory
    if not Path(outdir).exists():
        Path(outdir).mkdir()

    # convert counts to R matrix
    counts_flat = []
    for l in counts.transpose().to_numpy().tolist():
        for i in l:
            counts_flat.append(i)
    counts_mat = ro.r.matrix(counts_flat, nrow=counts.shape[0], ncol=counts.shape[1])

    # run polyester
    polyester = importr("polyester")
    polyester.simulate_experiment_countmat(
        fasta=txome_fa,
        readmat=counts_mat,
        outdir=outdir,
        paired=True,
        error_model="illumina5",
        gzip=True,
        seed=12,
    )


def run_polyester(txome_fa: str, counts: pd.DataFrame, outdir: str, n_jobs: int):
    """
    Run polyester to simulate reads from a transcriptome in parallel
    :param txome_fa: path to transcriptome fasta
    :param counts: pandas dataframe with reads per transcript per sample
    :param outdir: output directory
    :param n_jobs: number of parallel processes to use
    """
    # check input
    assert Path(txome_fa).exists(), f"No such file {txome_fa}"

    # make output directory
    if not Path(outdir).exists():
        Path(outdir).mkdir()

    # get full txome in memory
    txome = [tx for tx in SeqIO.parse(txome_fa, "fasta") if tx.id in counts.index]

    def producer(txome: list, counts: pd.DataFrame, chunk_size: int):
        """Split the counts and txome into chunks and return"""

        for i in range(0, counts.shape[0], chunk_size):
            this_counts = counts.iloc[i : i + chunk_size]
            # create a temporary fasta file with the transcripts in this chunk
            this_txome = [tx for tx in txome if tx.id in this_counts.index]
            assert (
                len(this_txome) == this_counts.shape[0]
            ), "Not all transcripts in counts were found in txome"
            SeqIO.write(this_txome, f"{outdir}/polyester_{i}.fa", "fasta")
            yield f"{outdir}/polyester_{i}.fa", this_counts, f"{outdir}/polyester_{i}"

    # run polyester in parallel
    chunk_size = ceil(counts.shape[0] / n_jobs)
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(_polyester)(fa, cts, odir)
        for fa, cts, odir in producer(txome, counts, chunk_size)
    )

    # merge the simulated reads
    logging.info("Merging simulated reads")
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
