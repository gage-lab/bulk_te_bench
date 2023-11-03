rule simulate_reads:
    input:
        txome_fa=rules.make_txome.output.fa,
        counts=lambda wc: config["txomes"][wc.txome]["simulations"][wc.sim]["counts"],
    output:
        reads=directory("results/{txome}/{sim}/reads"),
        counts="results/{txome}/{sim}/true_counts.tsv",
    conda:
        "simulate_reads.yaml"
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
