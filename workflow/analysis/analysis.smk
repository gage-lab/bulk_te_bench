def get_sample():
    return glob_wildcards(rules.salmon_quant_reads.output.quant_tx).sample


rule analyze:
    input:
        truth=rules.simulate_te_counts.output.counts,
        new_counts=expand(rules.salmon_quant_reads.output.quant_tx, sample=get_sample()),
    output:
        "results/{txome}/{tx_sim}/{te_sim}/txome_analysis.ipynb",
    conda:
        "analysis.yaml"
    log:
        notebook="results/{txome}/{tx_sim}/{te_sim}/txome_analysis.ipynb",
    notebook:
        "txome_analysis.py.ipynb"
