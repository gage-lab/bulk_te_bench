def get_strandedness(wc):
    if config["txomes"][wc.txome]["simulations"][wc.sim]["strand_specific"]:
        return "forward"
    else:
        return "no"


rule tetranscripts:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        genes_gtf=rules.make_txome.output.genes_gtf,
        rmsk_gtf=rules.make_txome.output.rmsk_gtf,
    output:
        tecount="results/{txome}/{tx_sim}/{te_sim}/tetranscripts/{sample}/TEtranscripts_out.cntTable",
        telocal="results/{txome}/{tx_sim}/{te_sim}/tetranscripts/{sample}/TElocal_out.cntTable",
    conda:
        "tetranscripts.yaml"
    log:
        tecount="results/{txome}/{tx_sim}/{te_sim}/tetranscripts/{sample}/TEtranscripts.log",
        telocal="results/{txome}/{tx_sim}/{te_sim}/tetranscripts/{sample}/TElocal.log",
    params:
        strandedness=get_strandedness,
        mode="multi",  # can be multi or uniq
    script:
        "tetranscripts.py"
