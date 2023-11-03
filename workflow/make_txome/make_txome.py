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

import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import pyranges as pr
from myutils.rmsk import read_rmsk
from pyroe import make_spliceu_txome
from snakemake.shell import shell

# unzip gtf and genome if necessary
for f in [snakemake.input.gencode_gtf, snakemake.input.genome_fa]:
    if ".gz" in f:
        out = f.replace(".gz", "")
        shell(f"gzip -dcf {f} > {out}")
snakemake.input.genome_fa = snakemake.input.genome_fa.replace(".gz", "")
snakemake.input.gencode_gtf = snakemake.input.gencode_gtf.replace(".gz", "")


# check if snakemake.params.chrom is a list
chrs = snakemake.params.chrs
if not isinstance(chrs, list):
    raise ValueError(chrs + " must be a list, fix config.yaml")

# parse and filter rmsk
rmsk_query = snakemake.params.rmsk_query
logger.info(
    f"Parsing rmsk, filtering for {rmsk_query} from chromosome(s) {' '.join(chrs)}"
)
rmsk = (
    read_rmsk(snakemake.input.rmsk_out)
    .query(rmsk_query)[["genoName", "genoStart", "genoEnd", "strand"]]
    .query("genoName in @chrs")
    .rename(
        columns={
            "genoName": "Chromosome",
            "genoStart": "Start",
            "genoEnd": "End",
            "strand": "Strand",
        }
    )
)

logger.info(f"Extracting sequences from {snakemake.input.genome_fa}")
# open temp files for txome creation
rmsk_fa = NamedTemporaryFile(suffix=".fa")
tmp_fa = NamedTemporaryFile(suffix=".fa")

# use bedtools to extract sequences from fasta, save to new fasta
logger.info(f"Extracting {rmsk_query} sequences from {snakemake.input.genome_fa}")
pr.PyRanges(rmsk).to_bed(snakemake.output.rmsk)
shell(
    f"bedtools getfasta -s -fi {snakemake.input.genome_fa} -bed {snakemake.output.rmsk} -fo {rmsk_fa.name}"
)

# make fasta for this chromosome
logger.info(
    f"Extracting chromosome(s) {' '.join(chrs)} from {snakemake.input.genome_fa}"
)
shell(f"samtools faidx {snakemake.input.genome_fa}")
shell(f"samtools faidx {snakemake.input.genome_fa} {' '.join(chrs)} > {tmp_fa.name}")

# make gtf for this chromosome
logger.info(
    f"Extracting transcripts from {snakemake.input.gencode_gtf} from chromosome(s) {' '.join(chrs)}"
)

# TODO: put this in config??
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


gtf = (
    pr.read_gtf(
        snakemake.input.gencode_gtf, as_df=True, duplicate_attr=True
    )  # duplicate_attr=True to keep all tags for each record
    .query("Chromosome in @chrs")
    .loc[filter_gtf]
)
pr.PyRanges(gtf).to_gtf(snakemake.output.genes_gtf)

# make rmsk 1 based
rmsk.Start = rmsk.Start - 1
joint = pd.concat([gtf, rmsk])
pr.PyRanges(joint).to_gtf(snakemake.output.joint_gtf)

# make splicu transcriptome
logger.info(f"Making spliceu transcriptome for chromosome(s) {' '.join(chrs)}")
make_spliceu_txome(
    genome_path=tmp_fa.name,
    gtf_path=snakemake.output.genes_gtf,
    output_dir=str(Path(snakemake.output.fa).parent),
    filename_prefix="txome",
    bt_path=shutil.which("bedtools"),
    extra_spliced=rmsk_fa.name,
    dedup_seqs=True,
)

# close all temp files
rmsk_fa.close()
tmp_fa.close()
