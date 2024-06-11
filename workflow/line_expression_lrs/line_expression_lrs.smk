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
                rules.L1_detection.output,  # TODO: update this with last rule in this part of pipeline
                zip,
                sample=ss["sample"],
                libtype=ss.libtype,
            )


rule lrs_done:
    input:
        get_lrs_output,
