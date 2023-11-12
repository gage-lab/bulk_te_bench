def get_estimates(wc):
    "Collect all estimates for a given quantifier and dataset"
    if wc.sim == "real_ont":
        ss = pd.read_csv(config["txomes"][wc.txome]["ont_samplesheet"], sep="\t")
        return expand(
            longread_quantifiers[wc.quant],
            zip,
            sample=ss["sample"],
            libtype=ss["libtype"],
            replicate=ss["replicate"],
            allow_missing=True,
        )
    else:
        my_wc = {}
        my_wc["tx_sim"], my_wc["te_sim"] = wc.sim.split("/")
        my_wc["txome"] = wc.txome
        checkpt_output = checkpoints.simulate_reads.get(**my_wc).output[0]
        samples = glob_wildcards(
            os.path.join(checkpt_output, "{sample}_1.fasta.gz")
        ).sample

        return expand(quantifiers[wc.quant], sample=samples, allow_missing=True)


rule aggregate:
    input:
        get_estimates,
    output:
        "results/{txome}/{sim}/{quant}_counts.tsv",
    conda:
        "aggregate.yaml"
    log:
        "results/{txome}/{sim}/{quant}_aggregate.log",
    script:
        "aggregate.py"
