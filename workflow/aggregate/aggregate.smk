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
    elif wc.sim == "real_illumina":
        ss = pd.read_csv(config["txomes"][wc.txome]["illumina_samplesheet"], sep="\t")
        return expand(
            shortread_quantifiers[wc.quant], sample=ss["sample"], allow_missing=True
        )
    else:
        my_wc = {
            "txome": wc.txome,
            "sim": wc.sim,
        }
        checkpt_output = checkpoints.concat_txte_simulations.get(**my_wc).output.reads
        samples = glob_wildcards(
            os.path.join(checkpt_output, "{sample}_1.fasta.gz")
        ).sample

        if len(samples) == 0:
            print(
                f"No samples found for {wc.txome} {wc.sim}: checkpt_output={checkpt_output}"
            )
            raise ValueError(f"No samples found for {wc.txome} {wc.sim}")

        return expand(
            shortread_quantifiers[wc.quant], sample=samples, allow_missing=True
        )


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
