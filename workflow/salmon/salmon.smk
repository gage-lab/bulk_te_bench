rule salmon_index:
    input:
        rules.make_txome.output.txome_fa,
    output:
        multiext(
            "results/{txome}/salmon_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
    log:
        "results/{txome}/salmon_index.log",
    threads: 2
    params:
        k=31,
    conda:
        "salmon.yaml"
    shell:
        "salmon index -t {input} -i $(dirname {output}) -k {params.k} > {log} 2>&1"


rule salmon_quant_reads:
    input:
        unpack(get_fq),
        gtf=rules.make_txome.output.joint_gtf,
        index=rules.salmon_index.output,
    output:
        quant_tx="results/{txome}/{sim}/salmon_quant_reads/{sample}/quant.sf",
        quant_ge="results/{txome}/{sim}/salmon_quant_reads/{sample}/quant.genes.sf",
        lib="results/{txome}/{sim}/salmon_quant_reads/{sample}/lib_format_counts.json",
        meta_info="results/{txome}/{sim}/salmon_quant_reads/{sample}/aux_info/meta_info.json",
        eq_classes="results/{txome}/{sim}/salmon_quant_reads/{sample}/aux_info/eq_classes.txt.gz",
    log:
        "results/{txome}/{sim}/salmon_quant_reads/{sample}/{sample}.log",
    params:
        libtype="A",
        extra="--gcBias --posBias --seqBias",
    threads: 2
    conda:
        "salmon.yaml"
    shell:
        """
        salmon quant \
            --threads {threads} \
            --geneMap {input.gtf} \
            --libType {params.libtype} \
            --dumpEq \
            --numGibbsSamples 30 \
            -i $(dirname {input.index[0]}) \
            -1 {input.fq1} \
            -2 {input.fq2} \
            -o $(dirname {output.quant_tx}) \
            {params.extra} > {log} 2>&1
        """


rule salmon_quant_bam:
    input:
        bam=rules.star_align.output.txome_bam,
        txome=rules.make_txome.output.txome_fa,
        gtf=rules.make_txome.output.joint_gtf,
    output:
        quant_tx="results/{txome}/{sim}/salmon_quant_bam/{sample}/quant.sf",
        quant_ge="results/{txome}/{sim}/salmon_quant_bam/{sample}/quant.genes.sf",
        meta_info="results/{txome}/{sim}/salmon_quant_bam/{sample}/aux_info/meta_info.json",
        eq_classes="results/{txome}/{sim}/salmon_quant_bam/{sample}/aux_info/eq_classes.txt.gz",
    log:
        "results/{txome}/{sim}/salmon_quant_bam/{sample}/{sample}.log",
    params:
        libtype="A",
        extra="--gcBias --posBias --seqBias",
    threads: 2
    conda:
        "salmon.yaml"
    shell:
        """
        salmon quant \
            --threads {threads} \
            --geneMap {input.gtf} \
            --libType {params.libtype} \
            --numGibbsSamples 30 \
            --dumpEq \
            -t {input.txome} \
            -a {input.bam} \
            -o $(dirname {output.quant_tx}) \
            {params.extra} > {log} 2>&1
        """


rule salmon_quant_bam_ont:
    input:
        bam=rules.minimap2.output,
        txome=rules.make_txome.output.txome_fa,
        gtf=rules.make_txome.output.joint_gtf,
    output:
        quant_tx="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}_rep{replicate}/quant.sf",
        quant_ge="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}_rep{replicate}/quant.genes.sf",
        meta_info="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}_rep{replicate}/aux_info/meta_info.json",
        eq_classes="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}_rep{replicate}/aux_info/eq_classes.txt.gz",
    log:
        "results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}_rep{replicate}/{sample}_{libtype}_rep{replicate}.log",
    params:
        libtype="A",
        extra="",  # cannot do bias correction with --ont
    threads: 2
    conda:
        "salmon.yaml"
    shell:
        """
        salmon quant \
            --threads {threads} \
            --geneMap {input.gtf} \
            --libType {params.libtype} \
            --numGibbsSamples 30 \
            --dumpEq \
            --ont \
            -t {input.txome} \
            -a {input.bam} \
            -o $(dirname {output.quant_tx}) \
            {params.extra} > {log} 2>&1
        """
