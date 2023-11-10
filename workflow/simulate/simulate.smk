def get_simulate_counts_input(wc):
    sim_config = config["txomes"][wc.txome]["simulations"][wc.tx_sim]
    out = {
        "txome_fa": rules.make_txome.output.txome_fa,
        "genes_gtf": rules.make_txome.output.genes_gtf,
        "rmsk_gtf": rules.make_txome.output.rmsk_gtf,
    }
    if wc.tx_sim == "gtex_sim":
        out["gtex_counts"] = sim_config["gtex_counts"]
        out["gtex_metadata"] = sim_config["gtex_metadata"]

    return out


rule simulate_tx_counts:
    input:
        unpack(get_simulate_counts_input),
    output:
        tx_counts="results/{txome}/{tx_sim}/tx_counts.tsv",
    log:
        "results/{txome}/{tx_sim}/simulate_tx_counts.log",
    conda:
        "simulate.yaml"
    script:
        "simulate_tx_counts.py"


rule simulate_te_counts:
    input:
        unpack(get_simulate_counts_input),
        tx_counts=rules.simulate_tx_counts.output.tx_counts,
    output:
        counts="results/{txome}/{tx_sim}/{te_sim}/true_counts.tsv",
        gene_counts="results/{txome}/{tx_sim}/{te_sim}/gene_counts.tsv",
    log:
        "results/{txome}/{tx_sim}/{te_sim}/simulate_te_counts.log",
    conda:
        "simulate.yaml"
    script:
        "simulate_te_counts.py"


checkpoint simulate_reads:
    input:
        txome_fa=rules.make_txome.output.txome_fa,
        counts=rules.simulate_te_counts.output.counts,
    output:
        reads=directory("results/{txome}/{tx_sim}/{te_sim}/reads"),
    conda:
        "simulate.yaml"
    threads: 32
    params:
        strand_specific=lambda wc: config["txomes"][wc.txome]["simulations"][
            wc.tx_sim
        ]["strand_specific"],
        readlen=lambda wc: config["txomes"][wc.txome]["simulations"][wc.tx_sim][
            "readlen"
        ],
    log:
        "results/{txome}/{tx_sim}/{te_sim}/simulate_reads.log",
    script:
        "simulate_reads.py"
