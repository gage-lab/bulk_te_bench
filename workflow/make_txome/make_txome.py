#!/usr/bin/env python
# Created on: Oct 24, 2023 at 9:17:15 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

from tempfile import NamedTemporaryFile

import pandas as pd
import pyranges as pr


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
    rmsktx["transcript_type"] = "retrogene"

    # exon level
    rmskex = rmsktx.copy().reset_index(drop=True)
    rmskex["Feature"] = "exon"
    rmskex["exon_id"] = rmskex.transcript_id
    rmskex["exon_number"] = 1

    return pd.concat([rmsk, rmsktx, rmskex]).sort_values(["Chromosome", "Start", "End"])


# from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build
BIOTYPES = [
    "protein_coding",
    "lncRNA",
    "lincRNA",
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


def make_unspliced_tx(gtf: pd.DataFrame):
    "Make unspliced transcripts from gtf, return gtf with unspliced transcripts added"

    # gene level
    ex = gtf[gtf.Feature == "exon"].set_index(["Chromosome", "Start", "End"]).copy()
    unsplicedtx = gtf[gtf.Feature == "gene"].copy()
    unsplicedtx.set_index(["Chromosome", "Start", "End"], inplace=True)
    n_unspliced = len(unsplicedtx)
    logger.info(f"Generating unspliced transcripts from {n_unspliced} genes")
    unsplicedtx = unsplicedtx[~unsplicedtx.index.isin(ex.index)].reset_index()
    n_unspliced -= len(unsplicedtx)
    logger.info(
        f"Removed {n_unspliced} genes that already have an unspliced transcript"
    )

    # tx level
    unsplicedtx["transcript_id"] = unsplicedtx.gene_id + "-I"
    unsplicedtx["transcript_type"] = unsplicedtx.gene_type + "_unspliced"
    unsplicedtx["Feature"] = "transcript"

    # exon level
    unsplicedex = unsplicedtx.copy()
    unsplicedex["exon_id"] = unsplicedex.gene_id + "-I"
    unsplicedex["exon_number"] = 1
    unsplicedex["Feature"] = "exon"

    # check for overlap with existing exons
    dup_ex = (
        gtf[gtf.Feature == "exon"][["Chromosome", "Start", "End"]]
        .isin(unsplicedex[["Chromosome", "Start", "End"]])
        .all(axis=1)
        .sum()
    )
    assert dup_ex == 0, f"{dup_ex} unspliced exons overlap with existing exons"

    return pd.concat([gtf, unsplicedtx, unsplicedex]).sort_values(
        ["Chromosome", "Start"]
    )


if __name__ == "__main__":

    from myutils.rmsk import read_rmsk
    from snakemake.shell import shell

    # send log to file
    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.INFO,
    )

    # unzip gtf and genome if necessary
    for f in [snakemake.input.gencode_gtf, snakemake.input.genome_fa]:
        if ".gz" in f:
            logger.info(f"Unzipping {f}")
            out = f.replace(".gz", "")
            shell(f"gzip -dcf {f} > {out}")
    snakemake.input.genome_fa = snakemake.input.genome_fa.replace(".gz", "")
    snakemake.input.gencode_gtf = snakemake.input.gencode_gtf.replace(".gz", "")

    # index and subset genome
    chrs = snakemake.params.chrs
    shell("samtools faidx {snakemake.input.genome_fa}")
    if len(chrs) == 0:
        # if not specified, get primary assembly chromosomes
        with open(snakemake.input.genome_fa + ".fai") as f:
            genome_chrs = [l.split("\t")[0] for l in f.readlines()]
        chrs = [c for c in genome_chrs if "_" not in c]
    my_chrs = " ".join(chrs)
    shell(
        "samtools faidx {snakemake.input.genome_fa} {my_chrs} > {snakemake.output.genome_fa}"
    )
    shell("samtools faidx {snakemake.output.genome_fa}")

    # parse and filter rmsk ###
    te_subfamilies = snakemake.params.te_subfamilies  # this is a list of subfams
    my_tes = " ".join(te_subfamilies)
    if len(te_subfamilies) == 0:
        my_tes = "Alu, SVA, and L1"

    logger.info(
        f"Parsing rmsk, filtering for full-length, transcriptionally-competent {my_tes} subfamilies on {my_chrs}"
    )
    rmsk = read_rmsk(snakemake.input.rmsk_out).query(
        "genoName in @chrs and has_promoter and is_full_length"
    )

    if len(te_subfamilies) > 0:
        rmsk = rmsk.query("repName in @te_subfamilies")

    rmsk.reset_index(drop=True, inplace=True)
    rmsk = rmsk_to_gtf(rmsk)

    ### parse and filter gencode ###
    logger.info(
        f"Parsing gencode GTF, filtering for high confidence transcripts on {my_chrs}"
    )

    genes = (
        pr.read_gtf(snakemake.input.gencode_gtf, as_df=True, duplicate_attr=True)
        .loc[filter_gtf]  # duplicate_attr=True to keep all tags for each record
        .query("Chromosome in @chrs")
        .reset_index(drop=True)
    )
    no_tx = genes.groupby("gene_id")["transcript_id"].count() == 0
    no_tx_genes = no_tx[no_tx].index
    genes = genes[genes["gene_id"].isin(no_tx_genes) == False]

    # get unspliced
    logger.info(f"Getting unspliced transcripts from {snakemake.input.gencode_gtf}")
    genes = make_unspliced_tx(genes).reset_index(drop=True)

    # annotate transcripts containing TEs and TEs contained in transcripts
    # use exons in genes dataframe to avoid transcripts that may not include the TE in a splice variant
    # Since unspliced transcripts are annotated with a single exon, if they contain a TE, the exon will contain it
    genes_with_tes = (
        pr.PyRanges(genes[genes.Feature == "exon"])
        .join(
            pr.PyRanges(rmsk[rmsk.Feature == "exon"]), report_overlap=True, suffix="_TE"
        )
        .df.query("Overlap == (End_TE - Start_TE)")[
            ["transcript_id", "transcript_id_TE", "Strand_TE"]
        ]
    )

    genes_with_tes["transcript_id_TE"] = (
        genes_with_tes["transcript_id_TE"].astype(str)
        + " ("
        + genes_with_tes["Strand_TE"].astype(str)
        + ")"
    )
    genes_with_tes = (
        genes_with_tes[["transcript_id", "transcript_id_TE"]]
        .groupby("transcript_id")
        .agg(lambda x: ",".join(x))
    )
    genes.loc[genes.Feature == "transcript", "contained_TEs"] = genes[
        genes.Feature == "transcript"
    ].join(genes_with_tes, how="left", on="transcript_id")["transcript_id_TE"]

    tes_in_genes = (
        pr.PyRanges(rmsk[rmsk.Feature == "exon"])
        .join(
            pr.PyRanges(genes[genes.Feature == "exon"]),
            how="left",
            report_overlap=True,
            suffix="_gene",
        )
        .df.query("Overlap == (End - Start)")[
            ["transcript_id", "transcript_id_gene", "Strand_gene"]
        ]
    )
    tes_in_genes["transcript_id_gene"] = (
        tes_in_genes["transcript_id_gene"].astype(str)
        + " ("
        + tes_in_genes["Strand_gene"].astype(str)
        + ")"
    )
    tes_in_genes = (
        tes_in_genes[["transcript_id", "transcript_id_gene"]]
        .groupby("transcript_id")
        .agg(lambda x: ",".join(x))
    )
    rmsk.loc[rmsk.Feature.isin(["exon", "transcript"]), "contained_in"] = rmsk[
        rmsk.Feature.isin(["exon", "transcript"])
    ].join(tes_in_genes, how="left", on="transcript_id")["transcript_id_gene"]

    if snakemake.wildcards.txome == "test_txome":  # type: ignore
        genes = genes.query("gene_id == 'ENSG00000100154.15'")

    # save genes_gtf and rmsk_gtf separately, then jointly.
    # use gffread to clean up the gene and joint gtf
    # -F = keep all GFF attributes (for non-exon features)
    #  --keep-exon-attrs : for -F option, do not attempt to reduce redundant exon/CDS attributes
    logger.info(f"Saving {len(genes[genes.Feature == 'transcript'])} transcripts")
    with NamedTemporaryFile(suffix=".gtf") as f:
        pr.PyRanges(genes).to_gtf(f.name)
        shell(
            "gffread -TFE --keep-genes {f.name} > {snakemake.output.genes_gtf} 2>> {snakemake.log}"
        )

    logger.info(f"Saving {len(rmsk[rmsk.Feature == 'exon'])} TEs")
    with NamedTemporaryFile(suffix=".gtf") as f:
        pr.PyRanges(rmsk[rmsk.Feature == "exon"]).to_gtf(f.name)
        shell(
            "gffread -TFE --keep-exon-attrs {f.name} | grep '\texon\t' > {snakemake.output.rmsk_gtf} 2>> {snakemake.log}"
        )

    logger.info(f"Concatenating rmsk and gencode")
    gtf = pd.concat([genes, rmsk]).sort_values(["Chromosome", "Start"])
    with NamedTemporaryFile(suffix=".gtf") as f:
        pr.PyRanges(gtf).to_gtf(f.name)
        shell(
            "gffread -TFE --keep-genes {f.name} > {snakemake.output.joint_gtf} 2>> {snakemake.log}"
        )

    # use gffread to extract unspliced, spliced, and TE sequences
    logger.info(f"Extracting sequences from {snakemake.input.genome_fa}")
    shell(
        "gffread "
        "-w {snakemake.output.txome_fa} "
        "-g {snakemake.output.genome_fa} "
        "{snakemake.output.joint_gtf} >> {snakemake.log} 2>&1"
    )
    shell("samtools faidx {snakemake.output.txome_fa}")

    logger.info("Done")
