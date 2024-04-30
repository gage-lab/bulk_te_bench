# OUTLINE
# - Read in data/samples.tsv
# - For every bam path in column:
# 1. rename path to replace beginning with "s3://sg-nex-data-blow5/"
# 2. download bam file to working directory using aws cli through command line
# 3. turn bam file into bigwig file using deeptools bamCoverage
# 4. run generate_heatmap.sh with bigwig file as input through command line
# 5. delete bam file from working directory

# Hardcoded variables - depend on data/illumina_samples.tsv
OUTDIR = "illumina_files/"
AWS_PREFIX = "https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/"
SAMPLES_PREFIX = AWS_PREFIX + "data/sequencing_data_illumina/bam/genome/"

import subprocess
from pathlib import Path

import pandas as pd

# read in samples.tsv
samples = pd.read_csv("data/illumina_samples.tsv", sep="\t")
total = len(samples["genome_bam_path"])

for i, path in enumerate(samples["genome_bam_path"]):
    # rename paths
    base_name = path.split("/")[-2]
    aws_path = path.replace(AWS_PREFIX, "s3://sg-nex-data/").replace(
        path.split("/")[-1], ""
    )
    local_bam_path = path.replace(SAMPLES_PREFIX, OUTDIR)
    local_bw_path = path.replace(SAMPLES_PREFIX, OUTDIR).replace(".bam", ".bw")
    print(f"***Downloading file {i+1}/{total}***")

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

    # turn bam file into bigwig file using deeptools bamCoverage
    print(f"***Turning bam into BigWig***")
    # TODO: check if option to make 0s empty
    result = subprocess.run(
        [
            "bamCoverage",
            "-b",
            local_bam_path,
            "-o",
            local_bw_path,
            "-p",
            "8",
            "--skipNAs",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert Path(local_bw_path).exists(), (
        "BigWig file not created. Error: " + result.stderr
    )

    # run generate_heatmap.sh with bigwig file as input through command line
    subprocess.run(["bash", "generate_heatmap.sh", local_bw_path])
    assert Path("heatmaps/" + base_name + ".bw.png").exists(), "Heatmap not created."

    print(f"***Heatmap generated {i+1}/{total}***")
    # subprocess.run(["rm", "-rf", local_bam_path])
