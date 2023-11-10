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

import pickle as pkl
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


read_fx = {
    "telocal": read_tetranscripts,
    "tecount": read_tetranscripts,
    "salmon_quant_reads_tx": read_salmon_quant,
    "salmon_quant_reads_ge": read_salmon_quant,
    "salmon_quant_bam_tx": read_salmon_quant,
    "salmon_quant_bam_ge": read_salmon_quant,
    "l1em": read_l1em,
}

# read in the files

res = {}

# get the samples
q = list(snakemake.input.keys())[0]
samples = [Path(f).parent.name for f in snakemake.input[q]]

for q in snakemake.input.keys():
    logger.info(f"Aggregating {q} files")

    # initialize dataframe
    out = read_fx[q](snakemake.input[q][0])
    res = pd.DataFrame(columns=samples, index=out.index)

    # get sample counts
    for f in snakemake.input[q]:
        sample = Path(f).parent.name
        res.loc[out.index, sample] = read_fx[q](f).loc[out.index, "counts"]

    # write out the files
    res.to_csv(snakemake.output[q], sep="\t", index=True, header=True)
