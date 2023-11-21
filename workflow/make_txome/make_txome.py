#!/usr/bin/env python
# Created on: Oct 24, 2023 at 9:17:15 PM
__author__ = "Michael Cuoco"

import logging

# send log to file
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import pyranges as pr
from myutils.rmsk import read_rmsk
from snakemake.shell import shell

# check if snakemake.params.chrom is a list
chrs = snakemake.params.chrs
if not isinstance(chrs, list):
    raise ValueError(chrs + " must be a list, fix config.yaml")

# unzip gtf and genome if necessary
for f in [snakemake.input.gencode_gtf, snakemake.input.genome_fa]:
    if ".gz" in f:
        out = f.replace(".gz", "")
        shell(f"gzip -dcf {f} > {out}")
snakemake.input.genome_fa = snakemake.input.genome_fa.replace(".gz", "")
snakemake.input.gencode_gtf = snakemake.input.gencode_gtf.replace(".gz", "")

# subset genome
my_chrs = " ".join(chrs)
shell("samtools faidx {snakemake.input.genome_fa}")
shell(
    "samtools faidx {snakemake.input.genome_fa} {my_chrs} > {snakemake.output.genome_fa}"
)

### parse and filter rmsk ###
te_subfamilies = snakemake.params.te_subfamilies  # this is a list of subfams
my_tes = " ".join(te_subfamilies)
if len(te_subfamilies) == 0:
    my_tes = "all TEs"


logger.info(f"Parsing rmsk, filtering for {my_tes} from chromosome(s) {my_chrs}")

rmsk = read_rmsk(snakemake.input.rmsk_out).query(
    "genoName in @chrs and has_promoter and is_full_length"
)

if len(te_subfamilies) > 0:
    rmsk = rmsk.query("repName in @te_subfamilies")

rmsk.reset_index(drop=True, inplace=True)


def rmsk_to_gtf(rmsk: pd.DataFrame) -> pd.DataFrame:
    "Add columns to make a repeatmasker gtf"

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
    rmsktx = rmsk.copy().reset_index()
    rmsktx["Feature"] = "transcript"
    rmsktx["transcript_id"] = (
        rmsktx["gene_id"] + "_dup" + rmsktx.groupby("gene_id").cumcount().astype(str)
    )

    # exon level
    rmskex = rmsktx.copy().reset_index()
    rmskex["Feature"] = "exon"
    rmskex["exon_id"] = rmskex.transcript_id
    rmskex["exon_number"] = 1

    return pd.concat([rmsk, rmsktx, rmskex]).sort_values(["Chromosome", "Start", "End"])


rmsk = rmsk_to_gtf(rmsk)

### parse and filter gencode ###
logger.info(
    f"Parsing gencode GTF, filtering for high confidence transcripts from chromosome(s) {my_chrs}"
)
# from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build
BIOTYPES = [
    "protein_coding",
    "lncRNA",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_LV_gene",
    "IG_V_gene",
    "IG_V_pseudogene",
    "IG_J_pseudogene",
    "IG_C_pseudogene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene",
    "TR_V_pseudogene",
    "TR_J_pseudogene",
]


def filter_gtf(gtf):
    """
    Filter gtf for high confidence transcripts
    1. gene_type in BIOTYPES and not on a Y chrom PAR region
    2. not annotated on an artifactual duplicate region of the genome assembly.
    3. transcript_type in BIOTYPES
    4. not a readthrough transcript
    """
    out = []
    for x in gtf.itertuples():
        tags = [x.tag] if type(x.tag) == float else x.tag.split(",")

        if (x.gene_type in BIOTYPES) and ("PAR" not in tags):
            if not hasattr(x, "artif_dupl") or (type(x.artif_dupl) == float):
                if x.Feature != "gene":
                    if (x.transcript_type in BIOTYPES) and (
                        "readthrough_transcript" not in tags
                    ):
                        out.append(True)
                        continue
                else:
                    out.append(True)
                    continue

        out.append(False)

    return out


genes = (
    pr.read_gtf(
        snakemake.input.gencode_gtf, as_df=True, duplicate_attr=True
    )  # duplicate_attr=True to keep all tags for each record
    .query("Chromosome in @chrs")
    .loc[filter_gtf]
)

# get unspliced
logger.info(f"Getting unspliced transcripts from {snakemake.input.gencode_gtf}")


def make_unspliced_tx(gtf: pd.DataFrame):
    "Make unspliced transcripts from gtf, return gtf with unspliced transcripts added"

    # gene level
    tx = (
        gtf[gtf.Feature == "transcript"]
        .set_index(["Chromosome", "Start", "End"])
        .copy()
    )
    unsplicedtx = gtf[gtf.Feature == "gene"].copy()
    unsplicedtx.set_index(["Chromosome", "Start", "End"], inplace=True)
    n_unspliced = len(unsplicedtx)
    logger.info(f"Generating unspliced transcripts from {n_unspliced} genes")
    unsplicedtx = unsplicedtx[~unsplicedtx.index.isin(tx.index)].reset_index()
    n_unspliced -= len(unsplicedtx)
    logger.info(
        f"Removed {n_unspliced} genes that already have an unspliced transcript"
    )

    # tx level
    unsplicedtx["transcript_id"] = unsplicedtx.gene_id + "-I"
    unsplicedtx["Feature"] = "transcript"

    # exon level
    unsplicedex = unsplicedtx.copy()
    unsplicedex["exon_id"] = unsplicedex.gene_id + "-I"
    unsplicedex["exon_number"] = 1
    unsplicedex["Feature"] = "exon"

    return pd.concat([gtf, unsplicedtx, unsplicedex]).sort_values(
        ["Chromosome", "Start", "End"]
    )


genes = make_unspliced_tx(genes)

# save genes_gtf and rmsk_gtf separately for tetranscripts
pr.PyRanges(genes).to_gtf(snakemake.output.genes_gtf)
pr.PyRanges(rmsk[rmsk.Feature == "exon"]).to_gtf(snakemake.output.rmsk_gtf)

# concat rmsk and gencode
logger.info(f"Concatenating rmsk and gencode")
gtf = pd.concat([genes, rmsk]).sort_values(["Chromosome", "Start", "End"])
pr.PyRanges(gtf).to_gtf(snakemake.output.joint_gtf)

# use gffread to extract unspliced, spliced, and TE sequences
shell(
    "gffread "
    "-w {snakemake.output.txome_fa} "
    "-g {snakemake.input.genome_fa} "
    "{snakemake.output.joint_gtf} >> {snakemake.log} 2>&1"
)

logger.info(f"Extracting sequences from {snakemake.input.genome_fa}")
