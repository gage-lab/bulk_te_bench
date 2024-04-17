#!/usr/bin/env python
# Created on: Apr 5, 2024 at 2:02:22 PM
__author__ = ["Michael Cuoco", "Joelle Faybishenko"]

import logging
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import perf_counter

import pysam

logging.basicConfig(level=logging.INFO)


def get_top_alns(inbam: str, output_dir: str):
    """
    Filter a BAM file to keep only the top scoring alignment for each read.
    :param infile: str: path to input BAM file
    """

    top_alns, top_score = {}, {}

    logging.info(f"Reading {inbam} and filtering top alignments...")
    start = perf_counter()
    with pysam.AlignmentFile(inbam, "rb") as bam:
        for aln in bam:
            # skip those that are bad reads
            if not aln.has_tag("AS"):
                continue
            id = aln.query_name + str(aln.is_read1)
            # if the read has already been seen, check if the current alignment is better
            if aln.query_name in top_alns:
                # if the current alignment is better, replace the old alignment and score
                if aln.get_tag("AS") > top_score[id]:
                    top_alns[id] = [aln]
                    top_score[id] = aln.get_tag("AS")
                # if the current alignment is the same as the best, add it to the list
                elif aln.get_tag("AS") == top_score[id]:
                    top_alns[id].append(aln)
            # if the read has not been seen, add it to the dictionaries
            else:
                top_alns[id] = [aln]
                top_score[id] = aln.get_tag("AS")
    logging.info(
        f"Finished filtering top alignments in {perf_counter() - start:.2f} seconds."
    )

    tempbam = NamedTemporaryFile(dir=output_dir, suffix=".bam")
    unique_reads = str(Path(output_dir) / Path(inbam).name).replace(
        ".bam", "_unique_reads.txt"
    )
    logging.info(
        f"Writing top alignments to {tempbam.name} and unique-mapping read IDs to {unique_reads}..."
    )
    start = perf_counter()
    with open(unique_reads, "w") as out_reads:
        with pysam.AlignmentFile(inbam, "rb") as bam:
            with pysam.AlignmentFile(tempbam.name, "wb", header=bam.header) as out_bam:
                for alns in top_alns.values():
                    # if there is only one alignment, write the read ID to the unique-mapping file
                    if len(alns) == 1:
                        out_reads.write(alns[0].query_name + "\n")
                    # write the top alignments for each read to the output BAM file
                    for a in alns:
                        out_bam.write(a)
    logging.info(
        f"Finished writing top alignments in {perf_counter() - start:.2f} seconds."
    )

    outbam = str(Path(output_dir) / Path(inbam).name).replace(".bam", "_top.bam")
    logging.info(f"Sorting {tempbam.name} to {outbam} and indexing...")
    pysam.sort("-o", outbam, tempbam.name)
    tempbam.close()
    pysam.index(outbam)

    logging.info(f"Finished indexing {outbam}.")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="Filter a BAM file to keep only the top scoring alignment for each read."
    )
    parser.add_argument("inbam", help="Path to input BAM file")
    parser.add_argument(
        "outdir", help="Optional: path to output directory", default=False
    )
    args = parser.parse_args()

    if not args.outdir:
        args.outdir = str(Path(args.inbam).parent)
    get_top_alns(args.inbam, args.outdir)
