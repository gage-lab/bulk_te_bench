# Benchmarking TE quantification approaches

[![Tests](https://github.com/gage-lab/bulk_te_bench/actions/workflows/main.yaml/badge.svg)](https://github.com/gage-lab/bulk_te_bench/actions/workflows/main.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.0-brightgreen.svg)](https://snakemake.github.io)

## Repository structure

```bash
.
├── README.md
├── conda.yaml # conda environment specification
├── analysis/ # analysis of pipeline outputs
├── workflow/ # snakemake pipeline for benchmarking
├── config.yaml # snakemake pipeline configuration
├── src/ # local python package
└── setup.py # install local python package
```

## The pipeline

Testing

```bash
# Generate test data:
cd .test/ngs-test-data
snakemake rnaseq/ref/txome.chr22.gtf rnaseq/ref/genome.chr22.fa scrnaseq_10x_v3/ref/rmsk_chr22.out -c1 --use-conda
cd ../..

snakemake all -c1 --use-conda --directory .test --show-failed-logs --rerun-incomplete
```

Running

```bash
snakemake all -c16 --use-conda --show-failed-logs --rerun-incomplete
```

Manually download to `./resources/GTEX`:

1. Transcript Counts - https://storage.cloud.google.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz
2. Gene Counts - https://storage.cloud.google.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
3. Sample Metadata - https://storage.cloud.google.com/adult-gtex/annotations/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

## Before committing

```bash
pre-commit install
```
