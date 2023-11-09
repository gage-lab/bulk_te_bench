def get_simulate_counts_input(wc):
    sim_config = config["txomes"][wc.txome]["simulations"][wc.sim]
    out = {
        "txome_fa": rules.make_txome.output.txome_fa,
        "genes_gtf": rules.make_txome.output.genes_gtf,
        "rmsk_gtf": rules.make_txome.output.rmsk_gtf,
    }
    if wc.sim == "gtex_sim":
        out["gtex_counts"] = sim_config["gtex_counts"]
        out["gtex_metadata"] = sim_config["gtex_metadata"]

    return out


rule simulate_counts:
    input:
        unpack(get_simulate_counts_input),
    output:
        counts="results/{txome}/{sim}/true_counts.tsv",
    log:
        "results/{txome}/{sim}/simulate_counts.log",
    conda:
        "simulate.yaml"
    script:
        "simulate_counts.py"


checkpoint simulate_reads:
    input:
        txome_fa=rules.make_txome.output.txome_fa,
        counts=rules.simulate_counts.output.counts,
    output:
        reads=directory("results/{txome}/{sim}/reads"),
    conda:
        "simulate.yaml"
    threads: 16
    params:
        strand_specific=lambda wc: config["txomes"][wc.txome]["simulations"][wc.sim][
            "strand_specific"
        ],
        readlen=lambda wc: config["txomes"][wc.txome]["simulations"][wc.sim]["readlen"],
    log:
        "results/{txome}/{sim}/simulate_reads.log",
    script:
        "simulate_reads.py"
