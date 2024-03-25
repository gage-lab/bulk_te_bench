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

"""

# # # DOWNLOAD BIGWIG FILES


# read in samples.tsv
samples = pd.read_csv("data/samples.tsv", sep="\t")
total = len(samples["bigwig_path"])
for i, path in enumerate(samples["bigwig_path"]):
	# rename paths
	aws_path = path.replace("https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/", "s3://sg-nex-data/")
	local_path = path.replace("https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/genome_browser_data/bigwig/", "longread_files/")
	print(f"***Downloading file {i+1}/{total}***")

	# download bigwig file to working directory using aws cli through command line
	result = subprocess.run(["aws", "s3", "cp", "--no-sign-request", aws_path, "longread_files/"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	assert Path(local_path).exists(), "File not downloaded. Error: " + result.stderr

	print("***\n")


# At 86:

Traceback (most recent call last):
  File "/netapp/LOG-G4/jfaybishenko/bulk_te_bench/analysis/download_long.py", line 20, in <module>
    aws_path = path.replace("https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/", "s3://sg-nex-data/")
AttributeError: 'float' object has no attribute 'replace'
"""

# # # DOWNLOAD BAM FILES


# Hardcoded variables - depend on data/illumina_samples.tsv
OUTDIR = "longread_files/"
AWS_PREFIX = "https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/"
SAMPLES_PREFIX = AWS_PREFIX + "data/sequencing_data_ont/bam/genome/"

import subprocess
from pathlib import Path

import pandas as pd

# read in samples.tsv
samples = pd.read_csv("data/samples.tsv", sep="\t")
total = len(samples["genome_bam_path"])

for i, path in enumerate(samples["genome_bam_path"]):
    # rename paths
    base_name = path.split("/")[-2]
    aws_path = path.replace(AWS_PREFIX, "s3://sg-nex-data/").replace(
        path.split("/")[-1], ""
    )
    local_bam_path = path.replace(SAMPLES_PREFIX, OUTDIR)
    print(f"***Downloading file {i+1}/{total}: {base_name} ***")

    # download bam + index files
    result = subprocess.run(
        [
            "aws",
            "s3",
            "cp",
            "--no-sign-request",
            aws_path,
            OUTDIR + base_name,
            "--recursive",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert Path(local_bam_path).exists(), "File not downloaded. Error: " + result.stderr
