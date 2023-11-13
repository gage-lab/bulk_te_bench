rule minimap2_index:
    input:
        target=rules.make_txome.output.txome_fa,
    output:
        "results/{txome}/minimap2_index/txome.mmi",
    log:
        "results/{txome}/minimap2_index/index.log",
    params:
        extra="",  # optional additional args
    threads: 3
    wrapper:
        "v2.11.1/bio/minimap2/index"


rule minimap2:
    input:
        unpack(get_fq),
        target=rules.minimap2_index.output,  # can be either genome index or genome fasta
    output:
        "results/{txome}/{sim}/{sample}_{libtype}_rep{replicate}.bam",
    log:
        "results/{txome}/{sim}/{sample}_{libtype}_rep{replicate}.log",
    params:
        extra=lambda wc: "-ax map-ont -k14 -uf"
        if wc.libtype == "directRNA"
        else "-ax map-ont",
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 3
    wrapper:
        "v2.11.1/bio/minimap2/aligner"
