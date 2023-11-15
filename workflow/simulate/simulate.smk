def get_simulate_counts_input(wc):
    out = {
        "txome_fa": rules.make_txome.output.txome_fa,
        "genes_gtf": rules.make_txome.output.genes_gtf,
        "rmsk_gtf": rules.make_txome.output.rmsk_gtf,
    }
    if wc.txte_sim == "gtex_sim":
        for i in ["gtex_counts", "gtex_metadata"]:
            if i not in config:
                raise ValueError(
                    f"{i} not specified in config.yaml, " "but gtex_sim requested"
                )
        out["gtex_counts"] = config["gtex_counts"]
        out["gtex_metadata"] = config["gtex_metadata"]

    return out


rule simulate_counts:
    input:
        unpack(get_simulate_counts_input),
    output:
        counts="results/{txome}/txte_simulations/{txte_sim}/simulated_counts_{nsamples}samples.tsv",
    log:
        "results/{txome}/txte_simulations/{txte_sim}/simulated_counts_{nsamples}samples.log",
    conda:
        "simulate.yaml"
    params:
        nsamples=lambda wc: wc.nsamples,
    script:
        "simulate_counts.py"


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
        counts=rules.simulate_counts.output.counts,
    output:
        reads=directory(
            "results/{txome}/txte_simulations/{txte_sim}/simulated_reads_{nsamples}samples_{readlen}bp_{strand}"
        ),
    log:
        "results/{txome}/txte_simulations/{txte_sim}/simulated_reads_{nsamples}samples_{readlen}bp_{strand}.log",
    conda:
        "simulate.yaml"
    threads: 32
    params:
        strand_specific=lambda wc: wc.strand,
        readlen=lambda wc: wc.readlen,
    script:
        "simulate_reads.py"


def get_sims(wc):
    """
    get reads and counts for each simuluation
    pass parameters for simulation from config and pass to te and tx simulation runs
    """

    tx_sim = config["txomes"][wc.txome]["simulations"][wc.sim]["tx"]
    te_sim = config["txomes"][wc.txome]["simulations"][wc.sim]["te"]

    other_wcs = {
        "readlen": config["txomes"][wc.txome]["simulations"][wc.sim]["readlen"],
        "strand": config["txomes"][wc.txome]["simulations"][wc.sim]["strand_specific"],
        "nsamples": config["txomes"][wc.txome]["simulations"][wc.sim]["nsamples"],
    }

    return {
        "tx_counts": expand(
            rules.simulate_counts.output.counts,
            txte_sim=tx_sim,
            **other_wcs,
            allow_missing=True,
        ),
        "tx_reads": expand(
            rules.simulate_reads.output.reads,
            txte_sim=tx_sim,
            **other_wcs,
            allow_missing=True,
        ),
        "te_counts": expand(
            rules.simulate_counts.output.counts,
            txte_sim=te_sim,
            **other_wcs,
            allow_missing=True,
        ),
        "te_reads": expand(
            rules.simulate_reads.output.reads,
            txte_sim=te_sim,
            **other_wcs,
            allow_missing=True,
        ),
    }


checkpoint concat_txte_simulations:
    input:
        unpack(get_sims),
    output:
        counts="results/{txome}/{sim}/simulated_counts.tsv",
        reads=directory("results/{txome}/{sim}/reads"),
    log:
        "results/{txome}/{sim}/concat_txte_simulations.log",
    conda:
        "simulate.yaml"
    script:
        "concat_txte_simulations.py"
