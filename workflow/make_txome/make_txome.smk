rule make_txome:
    input:
        genome_fa=remote_or_local(config["genome_fa"]),
        gencode_gtf=remote_or_local(config["gencode_gtf"]),
        rmsk_out=remote_or_local(config["rmsk_out"]),
    output:
        genome_fa="results/{txome}/genome.fa",
        txome_fa="results/{txome}/txome.fa",
        genes_gtf="results/{txome}/txome_genes.gtf",
        rmsk_gtf="results/{txome}/txome_rmsk.gtf",
        joint_gtf="results/{txome}/txome_joint.gtf",
    conda:
        "make_txome.yaml"
    log:
        "results/{txome}/make_txome.log",
    params:
        te_subfamilies=lambda wc: config["txomes"][wc.txome].get("te_subfamilies", []),
        chrs=lambda wc: config["txomes"][wc.txome].get("chrs", []),
    script:
        "make_txome.py"


rule kmer_similarity:
    input:
        gtf=rules.make_txome.output.joint_gtf,
        txome_fa=rules.make_txome.output.txome_fa,
    output:
        tsv="results/{txome}/te_kmer_sim.tsv",
        # png="results/{txome}/te_kmer_sim.png"
    conda:
        "make_txome.yaml"
    log:
        "results/{txome}/te_kmer_sim.log",
    script:
        "kmer_similarity.py"
