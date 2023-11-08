#!/usr/bin/env python
# Created on: Nov 2, 2023 at 1:17:47 PM
__author__ = "Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pickle
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

import pyranges as pr
from snakemake.shell import shell
from TElocal_Toolkit.TEindex import TEfeatures

with NamedTemporaryFile(suffix=".locInd") as tmp_loc:
    # run TEcount
    logger.info(f"Running TEcount on {snakemake.input.bam}")
    # create output directory
    outdir = Path(snakemake.output.tecount).parent
    Path(outdir).mkdir(parents=True, exist_ok=True)
    stem = str(outdir / Path(snakemake.output.tecount).stem)

    # run
    shell(
        "TEcount "
        "--BAM {snakemake.input.bam} "
        "--GTF {snakemake.input.genes_gtf} "
        "--TE {snakemake.input.rmsk_gtf} "
        "--stranded {snakemake.params.strandedness} "
        "--mode {snakemake.params.mode} "
        "--project {stem} "
        "--sortByPos "
        " --verbose 3 2> {snakemake.log.tecount}"
    )

    # run TElocal
    logger.info(f"Running TElocal on {snakemake.input.bam}")
    # create index
    ind = TEfeatures()
    ind.build(snakemake.input.rmsk_gtf)
    with open(tmp_loc.name, "wb") as f:
        pickle.dump(ind, f)

    # create output directory
    outdir = Path(snakemake.output.telocal).parent
    Path(outdir).mkdir(parents=True, exist_ok=True)
    stem = str(outdir / Path(snakemake.output.telocal).stem)

    # run
    shell(
        "TElocal "
        "--BAM {snakemake.input.bam} "
        "--GTF {snakemake.input.genes_gtf} "
        "--TE {tmp_loc.name} "
        "--stranded {snakemake.params.strandedness} "
        "--mode {snakemake.params.mode} "
        "--project {stem} "
        "--sortByPos "
        " --verbose 3 2> {snakemake.log.telocal}"
    )

sys.stderr.close()