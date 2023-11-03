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

import os
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

import pyranges as pr
from snakemake.shell import shell

rmsk_gtf = pr.read_gtf(snakemake.input.rmsk_gtf, as_df=True).set_index(
    ["Chromosome", "Start", "End"]
)
txome_rmsk = pr.read_bed(snakemake.input.txome_rmsk, as_df=True)
txome_rmsk.Start = txome_rmsk.Start - 1
txome_rmsk.set_index(["Chromosome", "Start", "End"], inplace=True)
rmsk = rmsk_gtf.loc[
    rmsk_gtf.index.isin(txome_rmsk.index.intersection(rmsk_gtf.index))
].reset_index()

outdir = Path(snakemake.output[0]).parent.__str__()
if not os.path.exists(outdir):
    os.makedirs(outdir)

# set variables for shell command
stem = outdir + "/" + Path(snakemake.output[0]).stem
cmd = "TElocal" if "telocal" in snakemake.rule else "TEcount"

with NamedTemporaryFile(suffix=".gtf") as rmsk_tmp:
    pr.PyRanges(rmsk).to_gtf(rmsk_tmp.name)
    shell(
        "{cmd} "
        "--BAM {snakemake.input.bam} "
        "--GTF {snakemake.input.txome_genes} "
        "--TE {rmsk_tmp.name} "
        "--stranded {snakemake.params.strandedness} "
        "--mode {snakemake.params.mode} "
        "--project {stem} "
        "--sortByPos "
        " --verbose 3 2> {snakemake.log}"
    )

sys.stderr.close()
