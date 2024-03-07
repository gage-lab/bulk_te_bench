#!/usr/bin/env python
# Created on: Nov 9, 2023 at 10:29:37 PM
__author__ = "Michael Cuoco"


import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG,
)

logger = logging.getLogger(__name__)

from pathlib import Path

import pandas as pd


# define functions for reading in files
def read_tetranscripts(f):
    """Reads a tetranscripts file into a pandas dataframe"""
    return pd.read_csv(f, sep="\t", index_col=0, names=["name", "counts"], header=0)


def read_salmon_quant(f):
    """Reads a salmon quant file into a pandas dataframe"""
    return pd.read_csv(
        f, sep="\t", index_col=0, usecols=[0, 4], names=["name", "counts"], header=0
    )


def read_l1em(f):
    """Reads a L1EM file into a pandas dataframe"""
    return pd.read_csv(
        f, sep="\t", index_col=0, usecols=[0, 1], names=["name", "counts"], header=0
    )


def read_oarfish_quant(f):
    """Reads a Oarfish file into a pandas dataframe"""
    return pd.read_csv(
        f, sep="\t", index_col=0, usecols=[0, 2], names=["name", "counts"], header=0
    )


read_fx = {
    "telocal": read_tetranscripts,
    "tecount": read_tetranscripts,
    "salmon_quant_bam_ont_tx": read_salmon_quant,
    "salmon_quant_bam_ont_ge": read_salmon_quant,
    "salmon_quant_reads_tx": read_salmon_quant,
    "salmon_quant_reads_ge": read_salmon_quant,
    "salmon_quant_bam_tx": read_salmon_quant,
    "salmon_quant_bam_ge": read_salmon_quant,
    "oarfish_quant_bam_ont": read_oarfish_quant,
    "l1em": read_l1em,
}

# get the samples
samples = [Path(f).parent.name for f in snakemake.input]

logger.info(f"Aggregating {snakemake.wildcards.quant} files")

# initialize dataframe
out = read_fx[snakemake.wildcards.quant](snakemake.input[0])
res = pd.DataFrame(columns=samples, index=out.index)

# get sample counts
for f in snakemake.input:
    sample = Path(f).parent.name
    res.loc[out.index, sample] = read_fx[snakemake.wildcards.quant](f).loc[
        out.index, "counts"
    ]

# write out the files
res.to_csv(snakemake.output[0], sep="\t", index=True, header=True)
