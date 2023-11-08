#!/usr/bin/env python
# Created on: Oct 24, 2023 at 4:55:26 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

from abc import abstractmethod
from pathlib import Path
from subprocess import Popen

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class Benchmark:
    """Base class for benchmarking transcriptome quantification tools"""

    def __init__(
        self, reads_dir: str | Path, txome_dir: str | Path, outdir: str | Path
    ) -> None:
        self.reads_dir = Path(reads_dir)
        self.txome_dir = Path(txome_dir)
        self.outdir = Path(outdir)
        self.benchmark = None

        # ensure these dirs exist
        for dir in [self.reads_dir, self.txome_dir]:
            if not dir.exists():
                raise FileNotFoundError(f"No such file or directory {dir}")

        # if outdir doesnt exist then make it
        if not self.outdir.exists():
            self.outdir.mkdir(exist_ok=True)

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

    def plot_difference(self) -> None:
        """
        Plot difference of estimated vs true counts in boxplots
        """

        if self.benchmark is None:
            raise ValueError("Must run read_counts() first!")

        ax = plt.gca()
        sns.boxplot(
            data=self.benchmark[~self.benchmark.tx_id.str.contains("chr")],
            x="sample",
            y="log2_abs_difference",
            ax=ax,
        )
        sns.stripplot(
            data=self.benchmark[self.benchmark.tx_id.str.contains("chr")],
            x="sample",
            y="log2_abs_difference",
            ax=ax,
            color="red",
        )
        ax.text(
            x=0.5,
            y=1.1,
            s="Difference of Estimated Counts from True Counts",
            fontsize=12,
            weight="bold",
            ha="center",
            va="bottom",
            transform=ax.transAxes,
        )
        ax.text(
            x=0.5,
            y=1.05,
            s="Red points = L1 transcripts",
            fontsize=8,
            alpha=0.75,
            ha="center",
            va="bottom",
            transform=ax.transAxes,
        )

    def plot_l1hs(self) -> None:
        """
        Plot estimated vs true counts for L1 transcripts
        """

        if self.benchmark is None:
            raise ValueError("Must run read_counts() first!")

        plot_df = self.benchmark[self.benchmark.tx_id.str.contains("chr")]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

        sns.lineplot(
            data=plot_df,
            x="count_true",
            y="count_estimated",
            hue="tx_id",
            palette="colorblind",
            ax=ax1,
        )

        sns.lineplot(
            data=plot_df,
            x="count_true",
            y="difference",
            hue="tx_id",
            palette="colorblind",
            ax=ax2,
        )

        # remove legends
        ax1.legend_.remove()

        # move ax2 legend outside
        handles, labels = ax2.get_legend_handles_labels()
        ax2.legend(
            handles=handles,
            labels=labels,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            frameon=False,
        )


class BenchmarkSalmon(Benchmark):
    """
    Wrapper for Salmon benchmarking
    1. Make transcriptome from genome fasta, GTF, and additioanl sequences (optional)
    2. Simulate reads from the transcriptome
    """

    def __init__(
        self, reads_dir: str | Path, txome_dir: str | Path, outdir: str | Path
    ) -> None:

        super().__init__(reads_dir, txome_dir, outdir)

    def index(self, k: int = 31, extra: str = "") -> None:
        """
        Index the transcriptome with salmon
        """

        cmd = f"salmon index -t {self.txome_dir}/txome.fa -i {self.outdir}/salmon_index -k {k} {extra} > {self.outdir}/salmon_index.log 2>&1"
        # self.logger.info(f"Indexing transcriptome with salmon: {cmd}")
        print(f"Indexing transcriptome with salmon: {cmd}")
        exit_code = Popen(cmd, shell=True).wait()
        if exit_code != 0:
            raise RuntimeError("Salmon index failed!")

    def quant(self, n_jobs: int = 8, extra: str = "") -> None:

        if not (self.outdir / "salmon_index").exists():
            raise FileNotFoundError("Must run index() first!")

        self.quant_dir = Path(str(self.reads_dir) + "_quant")
        self.quant_dir.mkdir(exist_ok=True)
        r1_reads = sorted(self.reads_dir.glob("*_1.fasta.gz"))
        r2_reads = sorted(self.reads_dir.glob("*_2.fasta.gz"))

        for r1, r2 in zip(r1_reads, r2_reads):
            sample = "_".join(r1.stem.split("_")[0:2])
            if sample != "_".join(r2.stem.split("_")[0:2]):
                raise ValueError("Reads are not paired properly!")

            cmd = f"salmon quant -g {self.txome_dir}/txome_t2g.tsv -i {self.outdir}/salmon_index -l A -1 {r1} -2 {r2} -o {self.quant_dir}/{sample} -p {n_jobs} {extra} > {self.quant_dir}/{sample}.log 2>&1"
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

        self.benchmark = benchmark
        return benchmark


# TODO: add classes BenchmarkTEtranscripts, BenchmarkL1EM, and BenchmarkSQuIRE
# use Benchmark as base class
# must have methods index(), quant(), and read_counts()

'''
class BenchmarkTEtranscripts(Benchmark):
    """
    Wrapper for TEtranscripts benchmarking
    """

    def __init__(
        self, reads_dir: str | Path, txome_dir: str | Path, outdir: str | Path
    ) -> None:

        super().__init__(reads_dir, txome_dir, outdir)

    def make_gtf(self) -> None:
        # read in self.txome_dir/extra.fa
        headers = []
        with open(f"{self.txome_dir}/extra.fa", "r") as infile:
            for line in infile:
                if line.startswith(">"):  # Check for header line
                    headers.append(line.strip())

        # write out headers to a file
        df_values = {
            "chr": [],
            "dot1": [],
            "start": [],
            "end": [],
            "dot2": [],
            "strand": [],
            "dot3": [],
            "name": [],
        }
        for header in headers:
            df_values["chr"].append(header.split(":")[0].split(">")[1])
            df_values["dot1"].append(".")
            df_values["start"].append(header.split(":")[1].split("-")[0])
            df_values["end"].append(header.split(":")[1].split("-")[1].split("(")[0])
            df_values["dot2"].append(".")
            df_values["strand"].append(
                header.split(":")[1].split("-")[1].split("(")[1].split(")")[0]
            )
            df_values["dot3"].append(".")
            df_values["name"].append('gene_id "L1HS"; gene_name "L1HS";')
        df = pd.DataFrame(df_values)
        df.to_csv(f"{self.txome_dir}/te.gtf", sep="\t", index=False, header=False)

    def index(self) -> None:
        """
        Indexing is making bam file for tetranscripts
        """
        if not (self.txome_dir / "te.gtf").exists():
            self.make_gtf()

    def quant(self) -> None:
        """
        Quantify the transcriptome with TEtranscripts
        """
        # check if extra.gtf exists
        if not (self.txome_dir / "extra.gtf").exists():
            raise FileNotFoundError("Must run index() first!")

        self.quant_dir = Path(str(self.reads_dir) + "_quant")
        self.quant_dir.mkdir(exist_ok=True)
        r1_reads = sorted(self.reads_dir.glob("*_1.fasta.gz"))
        r2_reads = sorted(self.reads_dir.glob("*_2.fasta.gz"))

        for r1, r2 in zip(r1_reads, r2_reads):
            sample = "_".join(r1.stem.split("_")[0:2])
            if sample != "_".join(r2.stem.split("_")[0:2]):
                raise ValueError("Reads are not paired properly!")

        cmd = f"TEcount -b OUTPUT OF INDEX --GTF {self.txome_dir}/clean_gtf.gtf --TE {self.txome_dir}/te.gtf --outdir {self.txome_dir}"
        logger.info(f"Running TEtranscripts quant: {cmd}")
        exit_code = Popen(cmd, shell=True).wait()
        if exit_code != 0:
            raise RuntimeError("TEtranscripts quant failed!")

    def read_counts(self) -> pd.DataFrame:
        """
        Get read counts from TEtranscripts output
        """

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
'''
