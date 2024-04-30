"""
Run bedtools intersect for each bam file and each L1 family.
"""

import subprocess

import numpy as np
import pandas as pd
import pyranges as pr
import pysam
from myutils import rmsk

# Use the best alignment bam file or not
BEST_ALIGNMENT = True
ILLUMINA = False  # If illumina is true then best should be false

"""
Calls bedtools intersect function.

@param file_path: path to the bam file
@param family: the family of the L1 element
@param best_alignment: whether to use the best alignment bam file or not
@return: the output file path

"""


def run_bedtools_intersect(file_path: str, family: str, best_alignment: bool):

    bed_files = [
        f"{family}.bed",
        f"full_intergenic_{family}.bed",
        f"full_intronic_{family}.bed",
        f"truncate_intergenic_{family}.bed",
        f"truncate_intronic_{family}.bed",
    ]

    prefix = "/".join(file_path.split("/")[:-1])
    for bed_file in bed_files:

        cmd = [
            "bedtools",
            "intersect",
            "-abam",
            file_path,
            "-b",
            "../data/" + bed_file,
            "-bed",
            "-wa",
            "-f",
            "0.5",
        ]
        if best_alignment:
            # Remove the last 4 characters from the bed file name and add "_best.bed"
            bed_file = bed_file[:-4] + "_best.bed"
        output_file = f"{prefix}/{bed_file}"

        with open(output_file, "w") as outfile:
            subprocess.run(cmd, stdout=outfile)
    return output_file


### ###

# Get the L1 families
L1_families = ["L1HS", "L1PA2", "L1PA3", "L1PA6"]

if BEST_ALIGNMENT:
    bam_files = [
        "../longread_files/SGNex_MCF7_directcDNA_replicate3_run3/SGNex_MCF7_directcDNA_replicate3_run3_bestAS.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate1_run2/SGNex_MCF7_directcDNA_replicate1_run2_bestAS.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate4_run2/SGNex_MCF7_directcDNA_replicate4_run2_bestAS.bam",
    ]

else:
    bam_files = [
        "../longread_files/SGNex_MCF7_directcDNA_replicate3_run3/SGNex_MCF7_directcDNA_replicate3_run3.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate1_run2/SGNex_MCF7_directcDNA_replicate1_run2.bam",
        "../longread_files/SGNex_MCF7_directcDNA_replicate4_run2/SGNex_MCF7_directcDNA_replicate4_run2.bam",
    ]


if ILLUMINA:
    bam_files = [
        "../illumina_files/SGNex_MCF7_Illumina_replicate2_run1/SGNex_MCF7_Illumina_replicate2_run1.bam",
        "../illumina_files/SGNex_MCF7_Illumina_replicate3_run1/SGNex_MCF7_Illumina_replicate3_run1.bam"
        "../illumina_files/SGNex_MCF7_Illumina_replicate4_run1/SGNex_MCF7_Illumina_replicate4_run1.bam",
    ]


# Run bedtools intersect for each bam file
for bam in bam_files:
    print(f"*****Processing {bam}*****")
    for i, family in enumerate(L1_families):
        print(f"{i}: processing {family}")
        last_file = run_bedtools_intersect(bam, family, BEST_ALIGNMENT)
