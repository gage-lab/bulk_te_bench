import numpy as np
import pandas as pd
import pyranges as pr
from myutils import rmsk

# # Creating bed files for LINE-1 Families

## ADAPTED FROM make_txome.py

"""
ADAPTED FROM make_txome.py

"""


def rmsk_to_gtf(rmsk: pd.DataFrame) -> pd.DataFrame:
    "convert RepeatMasker output file to gtf"

    rmsk.rename(
        columns={
            "genoName": "Chromosome",
            "genoStart": "Start",
            "genoEnd": "End",
            "strand": "Strand",
        },
        inplace=True,
    )

    # gene level
    rmsk.Start = rmsk.Start - 1  # make 0 based
    rmsk["gene_id"] = rmsk["repName"]
    rmsk["family_id"] = rmsk["repFamily"]
    rmsk["class_id"] = rmsk["repClass"]
    rmsk["gene_name"] = rmsk["repName"] + ":TE"
    rmsk = rmsk[
        [
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "gene_id",
            "family_id",
            "class_id",
            "gene_name",
        ]
    ].copy()
    rmsk["Source"] = "RepeatMasker"
    rmsk["Feature"] = "gene"
    rmsk["gene_type"] = "retrogene"

    # tx level
    rmsktx = rmsk.copy().reset_index(drop=True)
    rmsktx["Feature"] = "transcript"
    rmsktx["transcript_id"] = (
        rmsktx["gene_id"] + "_dup" + rmsktx.groupby("gene_id").cumcount().astype(str)
    )

    # exon level
    rmskex = rmsktx.copy().reset_index(drop=True)
    rmskex["Feature"] = "exon"
    rmskex["exon_id"] = rmskex.transcript_id
    rmskex["exon_number"] = 1

    return pd.concat([rmsk, rmsktx, rmskex]).sort_values(["Chromosome", "Start", "End"])


"""
Turns rmsk file into bed file.
@param rmsk_file: pd.DataFrame
@param genes: pd.DataFrame
@param family: str in file
@return bed: pd.DataFrame for testing
"""


def gtf_to_bed(rmsk_file, genes, family):
    "convert rmsk to bed"

    bed = rmsk_file[rmsk_file["repName"] == family]
    bed.to_csv(f"../data/{family}.bed", sep="\t", index=False, header=False)
    # get rid of chr in chromosome names
    bed["genoName"] = bed["genoName"].str.replace("chr", "")

    truncated_l1 = rmsk_to_gtf(bed[bed["is_full_length"] == False])
    full_l1 = rmsk_to_gtf(bed[bed["is_full_length"] == True])

    ### Full length intronic and intergenic regions
    full_l1["Length_TE"] = full_l1["End"] - full_l1["Start"]

    full_l1 = (
        pr.PyRanges(full_l1, int64=True)
        .join(
            pr.PyRanges(genes, int64=True),
            how="left",
            report_overlap=True,
            suffix="_gene",
        )
        .df
    )

    full_l1[(full_l1["Overlap"] == full_l1["Length_TE"])].to_csv(
        f"../data/full_intronic_{family}.bed", sep="\t", index=False, header=False
    )
    full_l1[(full_l1["Overlap"] != full_l1["Length_TE"])].to_csv(
        f"../data/full_intergenic_{family}.bed", sep="\t", index=False, header=False
    )

    ### Truncated intronic and intergenic regions
    truncated_l1["Length_TE"] = truncated_l1["End"] - truncated_l1["Start"]
    truncated_l1 = (
        pr.PyRanges(truncated_l1, int64=True)
        .join(
            pr.PyRanges(genes, int64=True),
            how="left",
            report_overlap=True,
            suffix="_gene",
        )
        .df
    )

    truncated_l1[(truncated_l1["Overlap"] == truncated_l1["Length_TE"])].to_csv(
        f"../data/truncate_intronic_{family}.bed", sep="\t", index=False, header=False
    )

    truncated_l1[(truncated_l1["Overlap"] != truncated_l1["Length_TE"])].to_csv(
        f"../data/truncate_intergenic_{family}.bed", sep="\t", index=False, header=False
    )

    return bed


# ## Read in gene functions

gtf = pr.read_gtf("../data/gencode.v26.basic.annotation.gtf.gz").df
genes = []
for gene_id, df in gtf.groupby("gene_id"):
    assert (
        df["Strand"].unique().shape[0] == 1
    ), f"This gene: {gene_id} is on two strands"

    row = {
        "Chromosome": df["Chromosome"].unique()[0],
        "Start": min(df["Start"]),
        "End": max(df["End"]),
        "Strand": df["Strand"].unique()[0],
        "gene_id": gene_id,
        "gene_name": df["gene_name"].unique()[0],
    }
    genes.append(row)

# additional filtering for future methods
genes = pd.DataFrame(genes)
genes["Length_Gene"] = genes["End"] - genes["Start"]
genes["Chromosome"] = genes["Chromosome"].str.replace(
    "chr", ""
)  # get rid of chr in chromosome names


rmsk_file = rmsk.read_rmsk("../data/hg38.fa.out.gz")
rmsk_file["genoName"].unique()
# keep only standard chromosomes
rmsk_file = rmsk_file[rmsk_file["genoName"].isin(genes["Chromosome"].unique())]


# ## Make bed files for LINE-1 families

L1_families = ["L1HS", "L1PA2", "L1PA3", "L1PA6"]

for family in L1_families:
    gtf_to_bed(rmsk_file, genes, family)
