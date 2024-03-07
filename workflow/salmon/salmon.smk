# create a decoy transcriptome for salmon
rule salmon_decoy:
    input:
        transcriptome=rules.make_txome.output.txome_fa,
        genome=rules.make_txome.output.genome_fa,
    output:
        gentrome="results/{txome}/gentrome.fa",
        decoys="results/{txome}/decoys.txt",
    threads: 2
    log:
        "results/{txome}/decoys.txt",
    wrapper:
        "v3.3.2/bio/salmon/decoys"


# create a salmon index with the decoy transcriptome
rule salmon_index:
    input:
        gentrome=rules.salmon_decoy.output.gentrome,
        decoys=rules.salmon_decoy.output.decoys,
    output:
        multiext(
            "results/{txome}/resources/salmon_index/",
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
        "results/{txome}/resources/salmon_index.log",
    threads: 2
    params:
        k=31,
    conda:
        "salmon.yaml"
    shell:
        "salmon index -t {input.gentrome} -i $(dirname {output}) --decoys {input.decoys} -k {params.k} > {log} 2>&1"


# quantify short reads with salmon (mapping-based)
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
        mappings="results/{txome}/{sim}/salmon_quant_reads/{sample}/mappings.bam",
    log:
        "results/{txome}/{sim}/salmon_quant_reads/{sample}/{sample}.log",
    params:
        libtype="A",
        extra="--validateMappings --seqBias --gcBias --posBias",
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
            {params.extra} --writeMappings 2> {log} | samtools view -bS > {output.mappings} > {log}
        """


# quantify short reads with salmon (alignment-based)
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


# quantify long reads with salmon (alignment-based)
rule salmon_quant_bam_ont:
    input:
        bam=expand(rules.minimap2.output, ome="txome", allowing=True),
        txome=rules.make_txome.output.txome_fa,
        gtf=rules.make_txome.output.joint_gtf,
    output:
        quant_tx="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}/quant.sf",
        quant_ge="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}/quant.genes.sf",
        meta_info="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}/aux_info/meta_info.json",
        eq_classes="results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}/aux_info/eq_classes.txt.gz",
    log:
        "results/{txome}/{sim}/salmon_quant_bam_ont/{sample}_{libtype}/{sample}_{libtype}.log",
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


# quantify long reads with oarfish (alignment-based)
rule oarfish_quant_bam_ont:
    input:
        bam=rules.minimap2.output,
        txome=rules.make_txome.output.txome_fa,
    output:
        quant="results/{txome}/{sim}/oarfish_quant_bam_ont/{sample}_{libtype}/{sample}_{libtype}.quant",
        meta_info="results/{txome}/{sim}/oarfish_quant_bam_ont/{sample}_{libtype}/{sample}_{libtype}.meta_info.json",
        infreps="results/{txome}/{sim}/oarfish_quant_bam_ont/{sample}_{libtype}/{sample}_{libtype}.infreps.pq",
    log:
        "results/{txome}/{sim}/oarfish_quant_bam_ont/{sample}_{libtype}/{sample}_{libtype}.log",
    params:
        libtype="A",
        extra="",  # cannot do bias correction with --ont
    threads: 2
    conda:
        "salmon.yaml"
    shell:
        """
        PREFIX="$(dirname {output.quant})/$(basename -s .quant {output.quant})"
        oarfish --alignments {input.bam} --threads {threads} --num-bootstraps 30 \
            --output $PREFIX > {log} 2>&1
        """
