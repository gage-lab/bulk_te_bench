# Benchmarking TE quantification approaches

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

## Testing the pipeline


```bash
# Generate test data:
cd .test/ngs-test-data
snakemake rnaseq/ref/txome.chr22.gtf rnaseq/ref/genome.chr22.fa rnaseq/ref/rmsk.chr22.gtf scrnaseq_10x_v3/ref/rmsk_chr22.out -c1 --use-conda

cd ../..
# test pipeline through make_txome and simulate_reads
snakemake results/test_txome/test_sim/reads -c1 --use-conda --directory .test --show-failed-logs --rerun-incomplete
```

## Before committing

```bash
pre-commit install
```
