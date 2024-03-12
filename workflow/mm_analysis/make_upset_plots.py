import logging
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from matplotlib import pyplot as plt
from upsetplot import UpSet, from_contents, from_memberships

"""
reads bedfiles for categories in path and return dataframe compatible with upsetplot
@param base_path: str path to bedfiles
@param categories: str list of categories to read
@return: dataframe compatible with upsetplot
"""


def read_data(
    bed_files, seperate_unique_multi=False, reads_path=None
) -> pd.DataFrame():

    if seperate_unique_multi:
        with open(reads_path, "r") as f:
            unique_reads = f.read().split()

    ids = {}
    # get list of read_ids for each sample
    for bed in bed_files:

        # if bed file is empty
        if Path(bed).stat().st_size == 0:
            continue

        if not seperate_unique_multi:
            ids[bed] = (
                pd.read_csv(bed, sep="\t", header=None)[3].drop_duplicates().to_list()
            )
        else:
            # iterate through reads and determine if unique or not
            ids[bed + "_unique"] = []
            ids[bed + "_multi"] = []
            for read in (
                pd.read_csv(bed, sep="\t", header=None)[3].drop_duplicates().to_list()
            ):
                if read in unique_reads:
                    ids[bed + "_unique"].append(read)
                else:
                    ids[bed + "_multi"].append(read)

    return from_contents(ids)


def make_a_plot(paths, seperate_unique_multi=False, reads_path=None):
    data = read_data(paths, seperate_unique_multi, reads_path)
    UpSet(data, subset_size="auto", show_counts=True, sort_categories_by="input").plot()
    # plt.suptitle(path)
    plt.show()  # change to  save figure
    plt.clf()


if __name__ == "__main__":

    # send log to file
    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)


# make_a_plot(snakemake.input.bed_files)  # type: ignore
