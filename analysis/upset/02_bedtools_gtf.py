#!/usr/bin/env python

"""
Run bedtools intersect for each bam file and each L1 family.
"""

__author__ = ["Joelle Faybishenko"]

import subprocess

import numpy as np
import pandas as pd
import pyranges as pr
from myutils import rmsk


def run_bedtools_intersect(file_path: str, family: str, outdir: str = None):
    """
    Calls bedtools intersect function.

    @param file_path: path to sorted bam file
    @param family: the family of the L1 element
    @return: the output file path

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

        cmd = [
            "bedtools",
            "intersect",
            "-abam",
            file_path,
            "-b",
            "../data/gtf_categories/" + cat + ".gtf",
            "-bed",
            "-wa",
            "-wb",
            "-f",
            "0.5",
        ]

        output_file = f"{prefix}/{cat}.bed"

        with open(output_file, "w") as outfile:
            subprocess.run(cmd, stdout=outfile)
    return output_file


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Intersect reads in BAM file with gtf files to find which reads intersect with L1 families. Output is bed files in same directory as input BAM file."
    )
    parser.add_argument("inbam", help="Path to input BAM file")
    parser.add_argument(
        "--outdir",
        help="Optional: path to output directory",
        default=None,
        required=False,
    )
    # add optional boolean arg for L1 family "chr"
    parser.add_argument(
        "--chr",
        help="Optional: add chr to L1 family name",
        action="store_true",
        default=False,
        required=False,
    )

    args = parser.parse_args()

    # Get the L1 families
    L1_families = ["L1HS", "L1PA2", "L1PA3", "L1PA6"]

    if args.chr:
        L1_families = [f"{family}chr" for family in L1_families]

    print(f"*****Processing {args.inbam}*****")
    for i, family in enumerate(L1_families):
        print(f"{i}: processing {family}")
        last_file = run_bedtools_intersect(args.inbam, family, args.outdir)
