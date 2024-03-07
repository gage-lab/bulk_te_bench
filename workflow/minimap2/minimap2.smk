def get_ome(wc):
    if wc.ome == "txome":
        return rules.make_txome.output.txome_fa
    else:
        return remote_or_local(config["genome_fa"])


rule minimap2_index:
    input:
        target=get_ome,
    output:
        "results/{txome}/resources/minimap2_index/{ome}.mmi",
    log:
        "results/{txome}/resources/minimap2_index/index.log",
    params:
        extra="",  # optional additional args
    threads: 3
    wrapper:
        "v2.11.1/bio/minimap2/index"


def get_minimap2_params(wildcards):
    if wildcards.libtype == "directRNA":
        return "-ax map-ont -uf -k14"
    else:
        return "-ax map-ont"


rule minimap2:
    input:
        unpack(get_fq),
        target=rules.minimap2_index.output,  # can be either genome index or genome fasta
    output:
        "results/{txome}/{sim}/{sample}_{libtype}_{ome}.bam",
    log:
        "results/{txome}/{sim}/{sample}_{libtype}_{ome}.log",
    params:
        extra=get_minimap2_params,
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 3
    wrapper:
        "v2.11.1/bio/minimap2/aligner"
