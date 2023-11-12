def get_strandedness(wc):
    tx_sim = wc.sim.split("/")[0]
    if config["txomes"][wc.txome]["simulations"][tx_sim]["strand_specific"]:
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
        tecount="results/{txome}/{sim}/tetranscripts/{sample}/TEtranscripts_out.cntTable",
        telocal="results/{txome}/{sim}/tetranscripts/{sample}/TElocal_out.cntTable",
    conda:
        "tetranscripts.yaml"
    log:
        tecount="results/{txome}/{sim}/tetranscripts/{sample}/TEtranscripts.log",
        telocal="results/{txome}/{sim}/tetranscripts/{sample}/TElocal.log",
    params:
        strandedness=get_strandedness,
        mode="multi",  # can be multi or uniq
    script:
        "tetranscripts.py"
