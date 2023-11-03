rule make_txome:
    input:
        genome_fa=remote_or_local(config["genome_fa"]),
        gencode_gtf=remote_or_local(config["gencode_gtf"]),
        rmsk_out=remote_or_local(config["rmsk_out"]),
    output:
        fa="results/{txome}/txome.fa",
        gencode_gtf="results/{txome}/txome_genes.gtf",
        rmsk="results/{txome}/txome_rmsk.bed",
    conda:
        "make_txome.yaml"
    log:
        "results/{txome}/make_txome.log",
    params:
        rmsk_query=lambda wc: config["txomes"][wc.txome]["rmsk_query"],
        chrs=lambda wc: config["txomes"][wc.txome]["chrs"],
    script:
        "make_txome.py"
