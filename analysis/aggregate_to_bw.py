import glob
import subprocess
from collections import defaultdict
from pathlib import Path

import pysam

OUTPUT_DIR = "longread_bigwigs/"


def parse_filename(filepath):
    # Extract the filename without the path and extension
    filename = filepath.split("/")[-1].split(".")[0]
    parts = filename.split("_")

    # Extract sample, replicate, and run information
    sample = "_".join(parts[:3])
    replicate = parts[-2]
    return sample, replicate


# Define the pattern to search for .bam files in the specified subdirectories
pattern = "longread_files/*/*.bam"

# Use glob.glob to find all files matching the pattern
bam_files = glob.glob(pattern)

# Group files by sample and replicate
grouped_files = defaultdict(list)
for file in bam_files:
    sample, replicate = parse_filename(file)
    grouped_files[(sample, replicate)].append(file)

# merge and turn into bigwigs
for key, files in grouped_files.items():
    output_bw = OUTPUT_DIR + key[0] + "_" + key[1] + ".bigwig"
    if Path(output_bw).exists():
        continue

    print(f"Processing {key[0]}_{key[1]}")
    if len(files) > 1:
        input_bam = files[0].replace("_run", "_")
        Path(input_bam).parents[0].mkdir(parents=True, exist_ok=True)
        pysam.merge(
            "-f", "-@", "8", "-o", input_bam, *files
        )  # Merge outputs a sorted bam
        pysam.index("-@", "8", input_bam)
    else:
        input_bam = files[0]

    result = subprocess.run(
        [
            "bamCoverage",
            "-b",
            input_bam,
            "-o",
            output_bw,
            "-p",
            "8",
            "--skipNAs",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert Path(output_bw).exists(), "BigWig file not created. Error: " + result.stderr
