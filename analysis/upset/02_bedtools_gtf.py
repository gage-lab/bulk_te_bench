"""
Run bedtools intersect for each bam file and each L1 family.
"""

import subprocess

import numpy as np
import pandas as pd
import pyranges as pr
import pysam
from myutils import rmsk

ILLUMINA = True


def run_bedtools_intersect(file_path: str, family: str):
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

    prefix = "/".join(file_path.split("/")[:-1])
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


### ###
# Get the L1 families
L1_families = ["L1HS", "L1PA2", "L1PA3", "L1PA6"]

"""bam_files = [
        "../longread_files/SGNex_MCF7_directcDNA_replicate3_run3/SGNex_MCF7_directcDNA_replicate3_run3.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate1_run2/SGNex_MCF7_directcDNA_replicate1_run2.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate4_run2/SGNex_MCF7_directcDNA_replicate4_run2.bam",
]


bam = "../longread_files/SGNex_MCF7_directcDNA_replicate3_run3/SGNex_MCF7_directcDNA_replicate3_run3.bam"
print(f"*****Processing {bam}*****")
for i, family in enumerate(L1_families):
    print(f"{i}: processing {family}")
    last_file = run_bedtools_intersect(bam, family, False)

"""

if ILLUMINA:
    bam_files = [
        "../illumina_files/SGNex_MCF7_Illumina_replicate2_run1/SGNex_MCF7_Illumina_replicate2_run1_top.bam",
        # "../illumina_files/SGNex_MCF7_Illumina_replicate3_run1/SGNex_MCF7_Illumina_replicate3_run1_top.bam",
        # "../illumina_files/SGNex_MCF7_Illumina_replicate4_run1/SGNex_MCF7_Illumina_replicate4_run1_top.bam",
    ]


for bam in bam_files:
    print(f"*****Processing {bam}*****")
    for i, family in enumerate(L1_families):
        print(f"{i}: processing {family}")
        last_file = run_bedtools_intersect(bam, family)
