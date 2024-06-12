rule download_line_expresssion_lrs:
    output:
        multiext(
            "results/LINE-Expression-LRS/scripts/",
            "01_preprocess_input.sh",
            "02_preprocess_mapping.sh",
            "03_L1_detection.sh",
            "04_map_hg38.sh",
            "05_map_qc_LRS.sh",
            "06_read_filter.sh",
            "07_exon_filter.sh",
            "08_L1_loci_filter.sh",
            "09_map_qc_LRS.sh",
            "10_normalization_wgt_avg.sh",
        ),
    shell:
        """
        git clone https://github.com/WGLab/LINE-Expression-LRS.git results/temp
        mv results/temp/* results/LINE-Expression-LRS
        rm -rf results/temp
        """


def get_lrs_fq(wc):
    for txome in config["txomes"]:
        if "ont_samplesheet" in config["txomes"][txome]:
            print(f"found ont samplesheet for {txome}")
            ss = pd.read_csv(config["txomes"][txome]["ont_samplesheet"], sep="\t")
            return ss.loc[
                (ss["libtype"] == wc.libtype) & (ss["sample"] == wc.sample), "fq"
            ].tolist()


rule preprocess_input:
    input:
        fastq=get_lrs_fq,
        script=rules.download_line_expresssion_lrs.output[0],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}.fasta",
    conda:
        "line_expression_lrs.yaml"
    params:
        libtype=lambda wc: wc.libtype,
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} <(gzip -dc ../../../{input.fastq}) FASTQ {params.libtype}
        """


rule preprocess_mapping:
    input:
        step1=rules.preprocess_input.output,
        script=rules.download_line_expresssion_lrs.output[1],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}_mapped_cDNA_1kb.sam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}_mapped_cDNA_1kb.fa",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample}
        """


rule L1_detection:
    input:
        step2=rules.preprocess_mapping.output[1],
        script=rules.download_line_expresssion_lrs.output[2],
    output:
        multiext(
            "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process/{sample}_{libtype}_mapped_cDNA_1kb.fa.",
            "align",
            "cat",
            "masked",
            "ori.out",
            "out",
            "out.xm",
            "tbl",
        ),
        "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process/{sample}_{libtype}_div10.fa",
        "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process/div10_LINEs.out",
        "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process/div10_readIDs.txt",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample}
        """


rule move_reference_genome:
    input:
        step1=rules.preprocess_input.output,
        ref=remote_or_local(config["genome_fa"]),
        gen=remote_or_local(config["gencode_gtf"]),
    output:
        newref="results/LINE-Expression-LRS/{sample}_{libtype}/references/hg38.fa",  # disclaimer: there is no guarantee these are hg38 and gencode v40
        newgen="results/LINE-Expression-LRS/{sample}_{libtype}/references/gencode.v40.annotation.bed",
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        mkdir -p results/LINE-Expression-LRS/{params.sample}/references
        cp {input.ref} {output.newref}
        cp {input.gen} {output.newgen}
        """


rule map_hg38:
    input:
        step3=rules.L1_detection.output[0],
        script=rules.download_line_expresssion_lrs.output[3],
        ref=rules.move_reference_genome.output[0],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample}
        """


rule map_qc_LRS:
    input:
        bam_input=rules.map_hg38.output[2],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}.log",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}/bam_summary.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}/img",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}/st_bam_statistics_dynamic.html",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}/st_bam_statistics.html",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:  # TODO add threads
        """
        cd results/LINE-Expression-LRS/{params.sample}/c_hg38_mapping_LRS/
        longreadsum bam -i $(basename {input.bam_input}) -o {params.sample}
        """


### after this we need a new parameter: active, inactive, or ORF2. Moving forwards with ORF2 for now for now
### I think this will not work because the OG authors had no idea how relative variables worked
rule read_filter:
    input:
        step5=rules.map_qc_LRS.output[1],
        script=rules.download_line_expresssion_lrs.output[5],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_read_filter_passed.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_read_filter_passed.sorted.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_read_filter_passed.sorted.bam.bai",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_bedgraph.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_bedgraph_clean.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/{l1_ref_type}/{sample}_{libtype}_bedgraph_sorted.bg",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        l1_ref_type="ORF2",  #l1_ref_type=lambda wc: wc.l1_ref_type
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} {params.l1_ref_type}
        """


def get_lrs_output(wc):
    for txome in config["txomes"]:
        if "ont_samplesheet" in config["txomes"][txome]:
            print(f"found ont samplesheet for {txome}")
            ss = pd.read_csv(config["txomes"][txome]["ont_samplesheet"], sep="\t")
            # NOTE: remove direct from libtype if present!
            ss["libtype"] = ss.libtype.apply(
                lambda x: x.lstrip("direct") if "direct" in x else x
            )
            return expand(
                rules.map_qc_LRS.output,  # TODO: update this with last rule in this part of pipeline
                zip,
                sample=ss["sample"],
                libtype=ss.libtype,
            )


rule lrs_done:
    input:
        get_lrs_output,
