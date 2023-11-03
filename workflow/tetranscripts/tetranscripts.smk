def get_strandedness(wc):
    if config["txomes"][wc.txome]["simulations"][wc.sim]["strand_specific"]:
        return "forward"
    else:
        return "no"


rule tetranscripts_count:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        txome_genes=rules.make_txome.output.genes_gtf,
        txome_rmsk=rules.make_txome.output.rmsk,
        rmsk_gtf=remote_or_local(config["rmsk_gtf"]),
    output:
        "results/{txome}/{sim}/tetranscripts/{sample}/TEtranscripts_out.cntTable",
    conda:
        "tetranscripts.yaml"
    shadow:
        "shallow"
    log:
        "results/{txome}/{sim}/tetranscripts/{sample}/TEtranscripts.log",
    params:
        strandedness=get_strandedness,
        mode="multi",  # can be multi or uniq
    script:
        "tetranscripts.py"
