#!/usr/bin/env bash
# Download resources for the project

cd "$(dirname "$0")"
mkdir -p resources

echo "Downloading hg38.fa"
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gzip -dc > resources/hg38.fa

echo "Download RepeatMasker annotation"
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz > resources/hg38.fa.out.gz

echo "Downloading transcript annotation"
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | gzip -dc > resources/gencode.v44.primary_assembly.basic.annotation.gtf
