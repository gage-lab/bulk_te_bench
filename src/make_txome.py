#!/usr/bin/env python
# Created on: Oct 24, 2023 at 1:32:22 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import requests, shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pyroe import make_spliceu_txome
from tempfile import NamedTemporaryFile
import pyranges as pr

REPBASE = "/iblm/logglun02/mcuoco/projects/salmonTE_testing/resources/hg38/humrep.ref"


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


def make_txome(
    outdir: str,
    genome_fa: str,
    txome_gtf: str,
    extra_fa: str = None,
    chromosome: str = None,
) -> None:
    """
    Make fasta file for transcriptome, including spliced and unsplied transcripts as well as TE sequences
    :param outdir: output directory
    :param genome_fa: fasta file of genome
    :param txome_gtf: transcriptome GTF file
    :param extra_fa: extra transcripts to add to transcriptome (optional)
    :param chromosome: chromosome to filter to (optional)
    """

    for f in [genome_fa, txome_gtf]:
        if not Path(f).exists():
            raise ValueError(f"No such file {f}")

    bt_path = shutil.which("bedtools")
    logger.info(f"Using bedtools from {bt_path}")

    with NamedTemporaryFile(suffix=".fa") as tmp_fa, NamedTemporaryFile(
        suffix=".gtf"
    ) as tmp_gtf:

        # make fasta for chromosome
        records = [s for s in SeqIO.parse(genome_fa, "fasta") if s.id == chromosome]
        SeqIO.write(records, tmp_fa.name, "fasta")
        logger.info(f"Chromosome {chromosome} fasta written to {tmp_fa.name}")

        # make gtf for chromosome
        gtf = pr.read_gtf(txome_gtf, as_df=True)
        pr.PyRanges(gtf[gtf["Chromosome"] == chromosome]).to_gtf(tmp_gtf.name)
        logger.info(f"Chromosome {chromosome} gtf written to {tmp_gtf.name}")

        logger.info(f"Saving spliceu txome to {outdir}")

        # extract isoforms and unspliced transcripts from genome
        # add TE sequences
        # save to FASTA
        make_spliceu_txome(
            tmp_fa.name,
            tmp_gtf.name,
            outdir,
            filename_prefix="txome",
            bt_path=bt_path,
            extra_spliced=extra_fa,
        )
