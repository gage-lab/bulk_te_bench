rule star_index:
    input:
        fasta=rules.make_txome.output.genome_fa,
    output:
        star_index=directory("results/{txome}/star_index"),
    threads: 8
    log:
        "results/{txome}/star_index.log",
    wrapper:
        "v1.20.0/bio/star/index"


rule star_align:
    input:
        fq1="results/{txome}/{sim}/reads/{sample}_1.fasta.gz",
        fq2="results/{txome}/{sim}/reads/{sample}_2.fasta.gz",
        idx=rules.star_index.output,
        gtf=rules.make_txome.output.joint_gtf,  # TODO: change this to gtf with rmsk and unspliced features
    output:
        genome_bam="results/{txome}/{sim}/star_align/{sample}/Aligned.out.bam",
        txome_bam="results/{txome}/{sim}/star_align/{sample}/Aligned.toTranscriptome.out.bam",
        log="results/{txome}/{sim}/star_align/{sample}/Log.out",
        log_final="results/{txome}/{sim}/star_align/{sample}/Log.final.out",
    threads: 8
    log:
        "results/{txome}/{sim}/star_align/{sample}/Log.err",
    params:
        # allowing for a maximum of 100 multi mapping loci and 200 anchors (used by Hammell Lab)
        # TODO: add description of each parameter
        extra=f"""--outSAMmultNmax -1 --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outMultimapperOrder Random --runRNGseed 777 --outSAMtype BAM Unsorted --sjdbScore 1""",
    conda:
        "star_align.yaml"
    script:
        # use custom script to return both genome_bam and txome_bam
        "star_align.py"


rule samtools_sort:
    input:
        rules.star_align.output.genome_bam,
    output:
        rules.star_align.output.genome_bam.replace("bam", "sorted.bam"),
    log:
        rules.star_align.log[0].replace("Log.err", "sort.log"),
    wrapper:
        "v1.20.0/bio/samtools/sort"


rule samtools_index:
    input:
        rules.samtools_sort.output,
    output:
        rules.samtools_sort.output[0] + ".bai",
    log:
        rules.samtools_sort.log[0].replace("sort", "index"),
    wrapper:
        "v1.20.0/bio/samtools/index"
