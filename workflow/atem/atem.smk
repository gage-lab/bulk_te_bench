rule atem:
    input:
        bam=rules.samtools_sort.output,
        bai=rules.samtools_index.output,
        genes_gtf=rules.make_txome.output.genes_gtf,
        rmsk_gtf=rules.make_txome.output.rmsk_gtf,
        meta_info=rules.salmon_quant_reads.output.meta_info,
    output:
        "results/{txome}/{sim}/atem/{sample}/out.txt",
    conda:
        "tetranscripts.yaml"
    log:
        "results/{txome}/{sim}/atem/{sample}/log.txt",
    script:
        "atem.py"
