import itertools
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam
from matplotlib import pyplot as plt
from upsetplot import UpSet, from_contents, from_memberships


def intersect(file_path: str, family: str, outdir: str = None):
    """
    Acts like runs bedtools intersect on a bam file with the gtf categories for an L1 family
    @param file_path: path to the bam file
    @param family: L1 family
    @return: True if successful
    """

    categories = [
        f"{family}",
        f"full_intergenic_{family}",
        f"full_intronic_{family}",
        f"truncate_intergenic_{family}",
        f"truncate_intronic_{family}",
    ]

    if outdir is None:
        prefix = "/".join(file_path.split("/")[:-1])
    else:
        prefix = outdir

    for cat in categories:
        output_file = f"{prefix}/{cat}.bed"
        gtf_file = "../data/gtf_categories/" + cat + ".gtf"

        # get list of transcript IDs from gtf file
        L1_loci = list(pr.read_gtf(gtf_file).df["transcript_id"])
        reads = []

        # iterate over bam file
        with pysam.AlignmentFile(file_path, "rb") as bam:
            for aln in bam:
                # get transctipt ID
                transcript_id = aln.reference_name
                # check if transcript ID is in list
                if transcript_id in L1_loci:
                    # write output of reads to bed file
                    reads.append(
                        [
                            aln.reference_name,
                            aln.reference_start,
                            aln.reference_end,
                            aln.query_name,
                            aln.is_reverse,
                        ]
                    )

        # write output of reads to bed file
        with open(output_file, "w") as f:
            for read in reads:
                f.write("\t".join(map(str, read)) + "\n")

    return True


def read_data(
    base_path: str,
    categories: list,
    seperate_unique_multi: bool = False,
    reads_path=None,
    illumina=False,
) -> pd.DataFrame():
    """
    reads bedfiles for categories in path and return dataframe compatible with upsetplot
    @param base_path: str path to bedfiles
    @param categories: str list of categories to read
    @return: dataframe compatible with upsetplot
    """

    if seperate_unique_multi:
        with open(reads_path, "r") as f:
            unique_reads = set(f.read().split())

    ids = {}
    # get list of read_ids for each sample
    for cat in categories:
        # if bed file is empty
        if Path(f"{base_path}/{cat}.bed").stat().st_size == 0:
            continue

        if illumina:
            readIDs = set(
                pd.read_csv(f"{base_path}/{cat}.bed", sep="\t", header=None)[3]
                .str[:-2]
                .drop_duplicates()
                .to_list()
            )
        else:
            readIDs = set(
                pd.read_csv(f"{base_path}/{cat}.bed", sep="\t", header=None)[3]
                .drop_duplicates()
                .to_list()
            )

        if not seperate_unique_multi:
            ids[cat] = readIDs

        else:
            # iterate through reads and determine if unique or not
            ids[cat + "_unique"] = []
            ids[cat + "_multi"] = []

            for read in readIDs:
                if read in unique_reads:
                    ids[cat + "_unique"].append(read)
                else:
                    ids[cat + "_multi"].append(read)

    return from_contents(ids)


def make_a_plot(
    path: str,
    categories: list,
    seperate_unique_multi: bool = False,
    reads_path=None,
    illumina=False,
):
    """
    uses upsetplot to make a plot of the data
    @param path: str path to bedfiles
    @param categories: str list of categories to read
    @param seperate_unique_multi: bool if True, will seperate unique and multi reads
    @param reads_path: str path to file with unique reads, only used if seperate_unique_multi is True
    @param illumina: bool if True, will remove last two characters from read ID
    """
    data = read_data(path, categories, seperate_unique_multi, reads_path, illumina)
    UpSet(data, subset_size="auto", show_counts=True, sort_categories_by="input").plot()
    plt.suptitle(path)
    plt.show()
    return plt
