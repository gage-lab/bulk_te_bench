# def get_intersect_bam_input(wc):


rule intersect_bam:
    output:
        bed=reads,
    conda:
        "mm_analysis.yaml"
    shell:
        """
        """
