rule make_txome:
    input:
        genome_fa=remote_or_local(config["genome_fa"]),
        gencode_gtf=remote_or_local(config["gencode_gtf"]),
        rmsk_out=remote_or_local(config["rmsk_out"]),
    output:
        genome_fa="results/{txome}/genome.fa",
        genome_fai="results/{txome}/genome.fa.fai",
        txome_fa="results/{txome}/txome.fa",
        txome_faidx="results/{txome}/txome.fa.fai",
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
        full_length=lambda wc: config["txomes"][wc.txome].get("full_length", False),
    script:
        "make_txome.py"


rule txome_report:
    input:
        gencode=remote_or_local(config["gencode_gtf"]),
        joint="results/{txome}/txome_joint.gtf",
    output:
        "results/{txome}/txome_report.ipynb",
    conda:
        "make_txome.yaml"
    log:
        notebook="results/{txome}/txome_report.ipynb",
    notebook:
        "txome_report.py.ipynb"


rule kmer_similarity:
    input:
        gtf=rules.make_txome.output.joint_gtf,
        txome_fa=rules.make_txome.output.txome_fa,
    output:
        npz="results/{txome}/te_kmer_sim.npz",
        # png="results/{txome}/te_kmer_sim.png"
    conda:
        "make_txome.yaml"
    log:
        "results/{txome}/te_kmer_sim.log",
    params:
        k=31,  # set this equal to the k used in salmon_index
    script:
        "kmer_similarity.py"
