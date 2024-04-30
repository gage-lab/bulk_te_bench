#!/bin/python3

# OUTLINE
# - Read in data/samples.tsv
# - For every bigwig path in column:
# 1. rename path to replace beginning with "s3://sg-nex-data-blow5/"
# 2. download bigwig file to working directory using aws cli through command line
# 3. run generate_heatmap.sh with bigwig file as input through command line
# 4. delete bigwig file from working directory

import subprocess
from pathlib import Path

import pandas as pd

# read in samples.tsv
samples = pd.read_csv("data/samples.tsv", sep="\t")
total = len(samples["bigwig_path"])
for i, path in enumerate(samples["bigwig_path"]):
    # rename paths
    aws_path = path.replace(
        "https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/", "s3://sg-nex-data/"
    )
    local_path = path.replace(
        "https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/genome_browser_data/bigwig/",
        "",
    )
    print(f"***Downloading file {i+1}/{total}***")

    # download bigwig file to working directory using aws cli through command line
    result = subprocess.run(
        ["aws", "s3", "cp", "--no-sign-request", aws_path, "."],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert Path(local_path).exists(), "File not downloaded. Error: " + result.stderr

    # run generate_heatmap.sh with bigwig file as input through command line
    subprocess.run(["bash", "generate_heatmap.sh", local_path])

    print(f"***Heatmap generated {i+1}/{total}***")
    subprocess.run(["rm", local_path])
    print("***\n***\n***\n")
