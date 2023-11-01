# Benchmarking TE quantification approaches

## Repository structure

```
.
├── README.md
├── conda.yaml # conda environment specification
├── download_resources.py # script to create resources directory with reference files
├── setup.py # install local python package
├── analysis # analysis-specific files
│   ├── l1hs_chr22.ipynb
│   ├── make_l1hs_chr22_txome.py
│   └── mikes_old_notebook.ipynb
└── src # local python package
    ├── __init__.py
    ├── txome.py
    └── bench.py
```

## Environment setup

```bash
# includes pip install -e . to install local python package from ./src/
mamba env create -f conda.yaml -p ./.conda

# download reference genome, transcript annotation, and parsed repeatmasker annotation
python download_resources.py
```

Manually download to `./resources/GTEX`:

1. Counts - https://storage.cloud.google.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz
2. Sample Metadata - https://storage.cloud.google.com/adult-gtex/annotations/v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

## Before committing

```bash
pre-commit install
```
