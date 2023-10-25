#!/usr/bin/env python
# Created on: Oct 24, 2023 at 1:32:22 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import shutil
from collections import defaultdict
from math import ceil
from pathlib import Path
from subprocess import Popen
from tempfile import NamedTemporaryFile

import pandas as pd
import pyranges as pr
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from pyroe import make_spliceu_txome
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

REPBASE = "/iblm/logglun02/mcuoco/projects/salmonTE_testing/resources/hg38/humrep.ref"

# TODO: move this as method of Txome?
def get_te_consensus(name: str, source: str, repbase: str | None = None) -> SeqRecord:
    """
    Get TE consensus sequence from Dfam or Repbase
    :param name: TE name
    :param source: "dfam" or "repbase"
    :param repbase: path to Repbase fasta file
    """

    if source == "dfam":
        # https://dfam.org/releases/Dfam_3.7/apidocs
        DFAM_URL = "https://dfam.org/api/families"
        logger.info(f"Retrieving {name} consensus sequence from DFAM API at {DFAM_URL}")

        params = {
            # The summary format is metadata-only and does not include
            "format": "full",
            "clade": "Homo sapiens",
            "name": name,
            "clade_relatives": "ancestors",
            "limit": 1,
        }

        r = requests.get(DFAM_URL, params=params).json()["results"][0]
        desc = f"{r['name']}\t{r['repeat_type_name']}\t{r['clades'][0].split(';')[-1]}"
        return SeqRecord(
            Seq(r["consensus_sequence"]), id=r["name"], name=r["name"], description=desc
        )

    elif source == "repbase":
        if not Path(repbase).exists():
            raise ValueError(f"no such file {repbase}")
        logger.info(f"Retrieving {name} consensus sequence from REPBASE at {repbase}")

        for s in SeqIO.parse(repbase, "fasta"):
            if s.id == name:
                return s
    else:
        raise ValueError("source must be dfam or repbase")


class Txome:
    """
    1. Make fasta file for transcriptome, including spliced and unsplied transcripts as well as TE sequences
    2. Index the transcriptome with salmon
    3. Cleanup transcript <-> gene map
    """

    def __init__(
        self,
        outdir: str | Path,
        genome_fa: str | Path,
        tx_gtf: str | Path,
        extra_fa: str | Path | None = None,
        chromosome: str = None,
    ) -> None:
        """
        :param outdir: output directory
        :param genome_fa: fasta file of genome
        :param txome_gtf: transcriptome GTF file
        :param extra_fa: extra transcripts to add to transcriptome (optional)
        :param chromosome: chromosome to filter to (optional)
        """
        # self.logger = logging.LoggerAdapter(logger, {"class": self.__class__.__name__})
        self.logger = logger
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.genome_fa = Path(genome_fa)
        self.extra_fa = Path(extra_fa) if extra_fa else None
        self.tx_gtf = Path(tx_gtf)
        # check inputs
        for f in [self.genome_fa, self.tx_gtf, self.extra_fa]:
            if f is not None and not f.exists():
                raise ValueError(f"No such file {f}")

        self.chromosome = chromosome

    def get_rmsk_seqs(self, parsed_rmsk: str | Path, query: str) -> None:
        """
        Get sequences from parsed rmsk file, save as fasta
        :param parsed_rmsk: path to parsed rmsk tsv file (by myutils.rmsk.read_rmsk())
        :param query: pandas query to filter on (e.g. "repName == 'L1HS' & length > 6000")
        """
        self.logger.info(f"Getting {query} sequences from {parsed_rmsk}")
        # NOTE: ignore has_promoter column for now, not sure if it is accurate
        rmsk = pd.read_csv(parsed_rmsk, sep="\t").query(query)
        rmsk = rmsk[["genoName", "genoStart", "genoEnd", "strand"]].rename(
            columns={
                "genoName": "Chromosome",
                "genoStart": "Start",
                "genoEnd": "End",
                "strand": "Strand",
            }
        )
        if self.chromosome:
            rmsk = rmsk[rmsk.Chromosome == self.chromosome]

        with NamedTemporaryFile(suffix=".bed") as bed:
            # save as bedfile
            pr.PyRanges(rmsk).to_bed(bed.name)
            # use bedtools to extract sequences from fasta, save to new fasta
            cmd = f"bedtools getfasta -s -fi {self.genome_fa} -bed {bed.name} -fo {self.outdir}/extra.fa"
            exit_code = Popen(cmd, shell=True).wait()
            if exit_code != 0:
                raise ValueError(f"bedtools getfasta failed with exit code {exit_code}")

        self.logger.info(f"Saved {len(rmsk)} sequences to {self.outdir}/extra.fa")
        if self.extra_fa:
            self.logger.warning(
                f"Extra fasta file {self.extra_fa} is being overwritten by {self.outdir}/extra.fa"
            )
        self.extra_fa = self.outdir / "extra.fa"

    def make_txome(self) -> None:

        bt_path = shutil.which("bedtools")
        self.logger.info(f"Using bedtools from {bt_path}")

        with NamedTemporaryFile(suffix=".fa") as tmp_fa, NamedTemporaryFile(
            suffix=".gtf"
        ) as tmp_gtf:

            # make fasta for chromosome
            records = [
                s
                for s in SeqIO.parse(self.genome_fa, "fasta")
                if s.id == self.chromosome
            ]
            SeqIO.write(records, tmp_fa.name, "fasta")
            self.logger.info(
                f"Chromosome {self.chromosome} fasta written to {tmp_fa.name}"
            )

            # make gtf for chromosome
            gtf = pr.read_gtf(self.tx_gtf, as_df=True)
            pr.PyRanges(gtf[gtf["Chromosome"] == self.chromosome]).to_gtf(tmp_gtf.name)
            self.logger.info(
                f"Chromosome {self.chromosome} gtf written to {tmp_gtf.name}"
            )

            # extract isoforms and unspliced transcripts from genome
            # add TE sequences
            # save to FASTA
            self.logger.info(f"Saving spliceu txome to {self.outdir}")
            make_spliceu_txome(
                tmp_fa.name,
                tmp_gtf.name,
                str(self.outdir),
                filename_prefix="txome",
                bt_path=bt_path,
                extra_spliced=self.extra_fa,
                dedup_seqs=True,
            )
            # reset root logger
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)

            self.txome_fa = self.outdir / "txome.fa"

    @staticmethod
    def _polyester_producer(
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

    def simulate_reads(self, counts: pd.DataFrame, name: str, n_jobs: int):
        """
        Run polyester to simulate reads from a transcriptome in parallel
        :param counts: pandas dataframe with reads per transcript per sample
        :param name: name of the run
        :param n_jobs: number of parallel threads to use
        """

        # make sure salmon index has been run
        if not (self.outdir / "salmon_index").exists():
            raise ValueError("salmon_index must be run before simulate_reads")

        outdir = self.outdir / name
        Path(outdir).mkdir(exist_ok=True)
        counts.to_csv(f"{outdir}/counts.tsv", sep="\t")

        n_transcripts, n_samples = counts.shape
        self.logger.info(
            f"Simulating reads from {n_transcripts} transcripts from {n_samples} samples with polyester"
        )

        # get full txome in memory
        txome = [
            tx for tx in SeqIO.parse(self.txome_fa, "fasta") if tx.id in counts.index
        ]

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
            for fa, cmat, odir in self._polyester_producer(
                txome, counts, chunk_size, outdir
            )
        )

        # merge the simulated reads
        self.logger.info("Merging simulated reads")
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
