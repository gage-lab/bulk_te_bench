#!/usr/bin/env python
# Created on: Oct 23, 2023 at 6:09:16 PM
__author__ = "Michael Cuoco"

from pathlib import Path
from subprocess import Popen
from myutils.rmsk import read_rmsk


def curl(url: str, outfile: str):
    "Download a file from a url"
    print(f"Downloading {url} to {outfile}")
    Popen(f"curl -s {url} | gzip -dc > {outfile}", shell=True).communicate()
    print("Done")
    Popen("sleep 1", shell=True).communicate()


if __name__ == "__main__":

    Path("resources").mkdir(exist_ok=True)
    curl(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "resources/hg38.fa",
    )
    curl(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz",
        "resources/gencode.v44.primary_assembly.basic.annotation.gtf",
    )
    read_rmsk(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz"
    ).to_csv("resources/hg38.rmsk.tsv", sep="\t", index=False)
