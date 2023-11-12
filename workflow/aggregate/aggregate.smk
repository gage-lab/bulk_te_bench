quantifiers = {
    "telocal": rules.tetranscripts.output.telocal,
    "tecount": rules.tetranscripts.output.tecount,
    "salmon_quant_reads_tx": rules.salmon_quant_reads.output.quant_tx,
    "salmon_quant_reads_ge": rules.salmon_quant_reads.output.quant_ge,
    "salmon_quant_bam_tx": rules.salmon_quant_bam.output.quant_tx,
    "salmon_quant_bam_ge": rules.salmon_quant_bam.output.quant_ge,
    "l1em": rules.l1em.output.full_counts,
}


def get_tx_te_sim(sim):
    tx_sim, te_sim = sim.split("/")
    return {"tx_sim": tx_sim, "te_sim": te_sim}


def get_estimates(wc):
    res = {}
    my_wc = get_tx_te_sim(wc.sim)
    my_wc["txome"] = wc.txome
    checkpt_output = checkpoints.simulate_reads.get(**my_wc).output[0]
    samples = glob_wildcards(os.path.join(checkpt_output, "{sample}_1.fasta.gz")).sample
    for n, q in quantifiers.items():
        res[n] = expand(q, sample=samples, allow_missing=True)
    return res


output = {}
for q in quantifiers.keys():
    output[q] = "results/{txome}/{sim}/" + q + "_counts.tsv"


rule aggregate:
    input:
        unpack(get_estimates),
    output:
        **output,
    conda:
        "aggregate.yaml"
    log:
        "results/{txome}/{sim}/aggregate.log",
    script:
        "aggregate.py"
