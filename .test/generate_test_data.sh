#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 11/11/23, 12:15 PM


# move in directory that this script is in
cd "$(dirname "$0")" && cd ngs-test-data

# check snakemake is installed
if ! command -v snakemake &> /dev/null
then
	echo "snakemake could not be found"
	exit
fi

snakemake rnaseq/ref/txome.chr22.gtf rnaseq/ref/genome.chr22.fa scrnaseq_10x_v3/ref/rmsk_chr22.out -c1 --use-conda

cd -
