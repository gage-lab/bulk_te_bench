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
bash .test/generate_test_data.sh

# test the pipeline
snakemake all -c1 --use-conda --directory .test --show-failed-logs
```

Running

```bash
# run the full pipeline
snakemake all -c16 --use-conda --show-failed-logs
```

Manually download to `./resources/GTEX`:

1. Gene Counts - https://storage.cloud.google.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
2. Sample Metadata - https://storage.cloud.google.com/adult-gtex/annotations/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

## Before committing

```bash
pre-commit install
```
