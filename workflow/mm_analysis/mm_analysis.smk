def get_intersect_bam_input(wc):
    out = {"txome": rules.make_txome.output.rmsk_gtf}

    if wc.best_AS:
        out["bam"] = rules.filter_bam.output.top_AS_bam
        out["unique_reads"] = rules.filter_bam.output.top_AS_unique_reads
    else:
        out["bam"] = ("results/{txome}/bam/{sample}/{sample}.bam",)
        out["unique_reads"] = rules.filter_bam.output.unique_reads


rule filter_bam:
    input:
        bam="results/{txome}/bam/{sample}/{sample}.bam",
    output:
        top_AS_bam="results/{txome}/bam/{sample}/{sample}_best_AS.bam",
        unique_reads="results/{txome}/bam/{sample}/unique_reads.txt",
        top_AS_unique_reads="results/{txome}/bam/{sample}/unique_reads_best_AS.txt.",
    conda:
        "mm_analysis.yaml"
    shell:
        """
        """


rule intersect_bam:
    input:
        unpack(get_intersect_bam_input),
        txome=rules.make_txome.output.rmsk_gtf,
        bam="results/{txome}/bam/{sample}/{sample}{best_AS}.bam",
    output:
        x="results/{txome}/salmon_quant_bam_ont/{sample}_{libtype}/quant.sf",
        bed=reads,
    conda:
        "mm_analysis.yaml"
    shell:
        """
        """


rule make_upset_plot:
    input:
        bed_files=rules.intersect_bam.output,
        reads="results/{txome}/bam/{sample}/unique_reads{best_AS}.txt",
    output:
        plot="results/{txome}/mm_analysis/{sample}/upset_plot{best_AS}.png",
    conda:
        "mm_analysis.yaml"
    script:
        "make_upset_plot.py"
