# TODO - log directory + symlink refs


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


rule move_reference_genome:
    input:
        step0=rules.download_line_expresssion_lrs.output,
        ref=remote_or_local(config["genome_fa"]),
        gen=remote_or_local(config["gencode_gtf"]),
    output:
        newref="results/LINE-Expression-LRS/references/hg38.fa",  # disclaimer: there is no guarantee these are hg38 and gencode v40
        newgen="results/LINE-Expression-LRS/references/gencode.v40.annotation.bed",
    shell:
        """
        cp {input.ref} {output.newref}
        cp {input.gen} {output.newgen}
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
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/preprocess_input.log",
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
    threads: 24
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/preprocess_mapping.log",
    shell:
        """
        logfile="../{params.sample}/log/preprocess_mapping.log"

        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} > $logfile 2>&1
        """


rule L1_detection:
    input:
        step2=rules.preprocess_mapping.output[1],
        script=rules.download_line_expresssion_lrs.output[2],
    output:
        multiext(
            "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process/{sample}_{libtype}_mapped_cDNA_1kb.fa.",
            "align",
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
    threads: 6
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/L1_detection.log",
    shell:
        """
        logfile="../{params.sample}/log/L1_detection.log"

        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} > $logfile 2>&1
        """


rule symbolic_link_refs:
    input:
        copied_ref=rules.move_reference_genome.output[0],
    output:
        new_refs=directory("results/LINE-Expression-LRS/{sample}_{libtype}/references/"),
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    shell:
        """
        full_source=$(realpath "results/LINE-Expression-LRS/references")
        full_target="$(realpath "results/LINE-Expression-LRS/{params.sample}")/references"

        ln -s "$full_source" "$full_target"
        """


rule map_hg38:
    input:
        step3=rules.L1_detection.output[0],
        script=rules.download_line_expresssion_lrs.output[3],
        ref=rules.symbolic_link_refs.output[0],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    threads: 8
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/map_hg38.log",
    shell:
        """
        logfile="../log/map_hg38.log"

        cd results/LINE-Expression-LRS/{params.sample}/c_hg38_mapping_LRS/

        FASTA_INPUT="../b_repeat_masker_process/{params.sample}_div10.fa"
        REF_SPLICE="../../references/gencode.v40.annotation.bed"
        REF_GENOME38="../../references/hg38.fa"


        echo "Mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome..." >> $logfile 2>&1

        minimap2 -ax splice --junc-bed $REF_SPLICE -uf --secondary=no -k14 -t 8 $REF_GENOME38 $FASTA_INPUT -o {params.sample}_hg38_mapped.sam >> $logfile 2>&1
        samtools view -Sb -o {params.sample}_hg38_mapped.bam {params.sample}_hg38_mapped.sam >> $logfile 2>&1
        samtools sort {params.sample}_hg38_mapped.bam -o {params.sample}_hg38_mapped.sorted_position.bam >> $logfile 2>&1
        samtools index {params.sample}_hg38_mapped.sorted_position.bam >> $logfile 2>&1

        echo "Finished mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome" >> $logfile 2>&1
        """


rule map_qc_LRS:
    input:
        bam_input=rules.map_hg38.output[2],
    output:
        directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}"
        ),
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/map_qc_LRS.log",
    shell:  # TODO add threads
        """
        logfile="../log/map_qc_LRS.log"

        cd results/LINE-Expression-LRS/{params.sample}/c_hg38_mapping_LRS/
        longreadsum bam -i $(basename {input.bam_input}) -o {params.sample} > $logfile 2>&1
        """


### after this we need a new parameter: active, inactive, or ORF2. Moving forwards with active for now for now
### I think this will not work because the OG authors had no idea how relative variables worked
rule read_filter:  # TODO turn every instance of "active" into a wc
    input:
        step5=rules.map_qc_LRS.output,
        script=rules.download_line_expresssion_lrs.output[5],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_read_filter_passed.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_read_filter_passed.sorted.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_read_filter_passed.sorted.bam.bai",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_bedgraph.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_bedgraph_clean.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/{sample}_{libtype}_bedgraph_sorted.bg",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        l1_ref_type="active",
    log:
        "../{sample}_{libtype}/log/active_map_qc_LRS.log",
    shell:
        """
        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} active > {log} 2>&1
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
