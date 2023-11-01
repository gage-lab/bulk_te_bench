#!/usr/bin/env python
# Created on: Oct 24, 2023 at 4:55:26 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

from abc import abstractmethod
from pathlib import Path
from subprocess import Popen

import numpy as np
import pandas as pd


class Benchmark:
    """Base class for benchmarking transcriptome quantification tools"""

    def __init__(self, reads_dir: str | Path, txome_dir: str | Path) -> None:
        self.reads_dir = Path(reads_dir)
        self.txome_dir = Path(txome_dir)

        # ensure these dirs exist
        for dir in [self.reads_dir, self.txome_dir]:
            if not dir.exists():
                raise FileNotFoundError(f"No such file or directory {dir}")

        pass

    @abstractmethod
    def index(self) -> None:
        """Index the transcriptome"""
        pass

    @abstractmethod
    def quant(self) -> None:
        """Quantify the transcriptome"""
        pass

    @abstractmethod
    def read_counts(self) -> pd.DataFrame:
        """Get read counts"""
        pass

    # TODO: add method plot_difference() to generate boxplot
    # TODO: add method plot_l1hs() to generate L1HS lineplots


class BenchmarkSalmon(Benchmark):
    """
    Wrapper for Salmon benchmarking
    1. Make transcriptome from genome fasta, GTF, and additioanl sequences (optional)
    2. Simulate reads from the transcriptome
    """

    def __init__(self, reads_dir: str | Path, txome_dir: str | Path) -> None:

        super().__init__(reads_dir, txome_dir)

    def index(self, k: int = 31, extra: str = "") -> None:
        """
        Index the transcriptome with salmon
        """

        cmd = f"salmon index -t {self.txome_dir}/txome.fa -i {self.outdir}/salmon_index -k {k} {extra} > {self.outdir}/salmon_index.log 2>&1"
        self.logger.info(f"Indexing transcriptome with salmon: {cmd}")
        exit_code = Popen(cmd, shell=True).wait()
        if exit_code != 0:
            raise RuntimeError("Salmon index failed!")

    def quant(self, n_jobs: int = 8, extra: str = "") -> None:

        if not (self.txome_dir / "salmon_index").exists():
            raise FileNotFoundError("Must run index() first!")

        self.quant_dir = Path(str(self.reads_dir) + "_quant")
        self.quant_dir.mkdir(exist_ok=True)
        r1_reads = sorted(self.reads_dir.glob("*_1.fasta.gz"))
        r2_reads = sorted(self.reads_dir.glob("*_2.fasta.gz"))

        for r1, r2 in zip(r1_reads, r2_reads):
            sample = "_".join(r1.stem.split("_")[0:2])
            if sample != "_".join(r2.stem.split("_")[0:2]):
                raise ValueError("Reads are not paired properly!")

            cmd = f"salmon quant -g {self.txome_dir}/txome_t2g.tsv -i {self.txome_dir}/salmon_index -l A -1 {r1} -2 {r2} -o {self.quant_dir}/{sample} -p {n_jobs} {extra} > {self.quant_dir}/{sample}.log 2>&1"
            logger.info(f"Running Salmon quant for {sample}")
            exit_code = Popen(cmd, shell=True).wait()
            if exit_code != 0:
                raise RuntimeError(f"Salmon quant failed for {sample}")

    def read_counts(self) -> pd.DataFrame:

        if not hasattr(self, "quant_dir"):
            raise ValueError("Must run quant() first!")

        # read in ground truth counts
        truth = pd.read_csv(self.reads_dir / "counts.tsv", sep="\t", index_col=0).melt(
            value_name="count", var_name="sample", ignore_index=False
        )
        truth["sample"] = truth["sample"].astype(int) + 1
        truth = truth.reset_index().set_index(["tx_id", "sample"])

        # read in estimated counts
        estimate = {}
        for f in self.quant_dir.rglob("quant.sf"):
            sample = int(f.parent.name.split("_")[1])
            estimate[sample] = pd.read_csv(f, sep="\t", index_col=0).iloc[:, 3]

        estimate = (
            pd.DataFrame(estimate)
            .melt(value_name="count", var_name="sample", ignore_index=False)
            .reset_index()
            .rename(columns={"Name": "tx_id"})
            .set_index(["tx_id", "sample"])
        )

        # combine truth and estimate, compute difference
        benchmark = truth.join(estimate, lsuffix="_true", rsuffix="_estimated")
        benchmark["difference"] = benchmark["count_true"] - benchmark["count_estimated"]
        benchmark["log2_abs_difference"] = np.log2(np.abs(benchmark["difference"]))
        benchmark.reset_index(inplace=True)

        return benchmark


# TODO: add classes BenchmarkTEtranscripts, BenchmarkL1EM, and BenchmarkSQuIRE
# use Benchmark as base class
# must have methods index(), quant(), and read_counts()
