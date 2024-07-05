
rule download_line_expresssion_lrs:
    output:
        multiext(
            "resources/LINE-Expression-LRS/scripts/",
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
        "resources/LINE-Expression-LRS/references/custom_LINE_reference.fasta",
    shell:
        """
        git clone https://github.com/WGLab/LINE-Expression-LRS.git resources/temp
        mv resources/temp/* resources/LINE-Expression-LRS
        rm -rf resources/temp
        """


# TODO : delete this rule and use ref/gen paths in future rules
rule move_reference_genome:
    input:
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
    output:
        dirA=directory("results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset"),
        dirB=directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/b_repeat_masker_process"
        ),
        dirC=directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS"
        ),
        dirD=directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification"
        ),
        fasta="results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}.fasta",
        fasta1kb="results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}_cDNA_1kb.fasta",
    conda:
        "line_expression_lrs.yaml"
    params:
        libtype=lambda wc: wc.libtype,
        filetype="FASTQ",
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/preprocess_input.log",
    shell:
        """
        echo "{params.sample}" >> {log}

        # Folder Preparation
        mkdir -p {output.dirA}
        mkdir -p {output.dirB}
        mkdir -p {output.dirC}
        mkdir -p {output.dirD}

        echo "All folders have been made!" >> {log}


        # Input File Preparation

        if [ {params.filetype} == "FASTQ" ]; then
            echo "Input file is:   FASTQ" >> {log}
            echo "Converting to FASTA..." >> {log}
            zcat {input.fastq} | awk 'NR%4==1{{printf ">%s\\n", substr($0,2)}} NR%4==2{{print}}' > {output.fasta} #TODO what if not gzipped


            echo "Filtering out reads less than 1kb..." >> {log}
            awk '/^>/ {{if (seqlen >= 1000) {{print header; print seq}} header=$0; seq=""; seqlen=0; next}} {{seq = seq $0; seqlen += length($0)}} END {{if (seqlen >= 1000) {{print header; print seq}}}}' {input.fastq} > {output.fasta1kb} # IDEALLY THIS OUTPUTS A TEMP/INTERM FILE WITH DIFF NAME

        fi

        # RNA to cDNA Conversion
        # If input file is RNA, convert to cDNA

        if [ {params.libtype} == "RNA" ]; then
            echo "Input file is:  RNA"  >> {log}
            echo " Converting to cDNA..."  >> {log}
            perl -pe 'tr/uU/tT/ unless />/' < {output.fasta1kb} > {output.fasta1kb} # HERE IS WHERE THE FILE IS ADJUSTED
        fi

        if [ {params.libtype} == "DNA" ]; then
            echo "Input file is:  DNA"  >> {log}
            #mv {output.fasta1kb} {output.fasta1kb} # OR HERE IS WHERE THE FILE IS RENAMED
        fi
        echo "Preprocessing complete!"  >> {log}
        echo""  >> {log}
        """


rule preprocess_mapping:
    input:
        fasta1kb=rules.preprocess_input.output.fasta1kb,
        ref_L1_mega="resources/LINE-Expression-LRS/references/custom_LINE_reference.fasta",
    output:
        sam="results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}_mapped_cDNA_1kb.sam",
        fa="results/LINE-Expression-LRS/{sample}_{libtype}/a_dataset/{sample}_{libtype}_mapped_cDNA_1kb.fa",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    threads: 24
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/preprocess_mapping.log",
    shell:
        """
        echo {params.sample} >> {log}

        echo "Mapping to Custom LINE Reference Library ..." >> {log}

        minimap2 -ax map-ont {input.ref_L1_mega} {input.fasta1kb} -t {threads} > {output.sam}
        samtools fasta {output.sam} -F 2308 > {output.fa} 2>{log}

        echo "Mapping to Custom LINE Reference Library complete!" >> {log}

        """


rule L1_detection:
    input:
        step2=rules.preprocess_mapping.output[1],
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
        #.cat/.cat.gz too
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


rule map_hg38:
    input:
        step3=rules.L1_detection.output[0],
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
rule read_filter:  # TODO turn every instance of "active" into a wc
    input:
        step5=rules.map_qc_LRS.output,
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph_clean.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph_sorted.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph.bg",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_L1_regions_reads.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/read_filter.log",
    shell:
        """
        logfile="../../../log/read_filter.log"
        sample_name={params.sample}
        L1_ref_type="active"

        echo "Read Filter on the $L1_ref_type Reference L1 Regions" > {log} 2>&1

        L1_ref_input="../../../../references/L1Base2_filtered/${{L1_ref_type}}_filtered.bed"
        sample_bam_input="../../../../${{sample_name}}/c_hg38_mapping_LRS/${{sample_name}}_hg38_mapped.sorted_position.bam"

        cd results/LINE-Expression-LRS/${{sample_name}}/d_LINE_quantification/

        mkdir -p $L1_ref_type; cd $L1_ref_type
        mkdir -p "read_filter"; cd "read_filter"

        # Output Files
        L1_regions_reads=${{sample_name}}_L1_regions_reads.bam
        echo "Generating filtered BAM file with reads only located within the L1 reference regions..." >> $logfile 2>&1
        samtools view -b -L "$L1_ref_input" -o "$L1_regions_reads" "$sample_bam_input"

        read_filter_bam=${{sample_name}}_read_filter_passed.bam
        echo "Removing reads with less than 90% of the read maps to the L1 reference regions" >> $logfile 2>&1
        bedtools intersect -a "$L1_regions_reads" -b "$L1_regions_reads" -f 0.9 > "$read_filter_bam"

        sorted_read_filter_bam=${{sample_name}}_read_filter_passed.sorted_position.bam
        echo "Sorting and Indexing the resulting Read Filter BAM file..." >> $logfile 2>&1
        samtools sort "$read_filter_bam" -o "$sorted_read_filter_bam"
        samtools index "$sorted_read_filter_bam"

        # Generate the bedgraph for the new BAM file that passed the Read Filter
        bedgraph_output=${{sample_name}}"_bedgraph.bg"
        echo "Generating the bedgraph..." >> $logfile 2>&1
        bedtools genomecov -ibam "$sorted_read_filter_bam" -bga -split > "$bedgraph_output"

        bedgraph_output_clean=${{sample_name}}"_bedgraph_clean.bg"
        echo "Cleaning the bedgraph..." >> $logfile 2>&1
        grep -v 'fix\|alt\|random\|[(]\|Un' $bedgraph_output > $bedgraph_output_clean

        bedgraph_sort_output=${{sample_name}}"_bedgraph_sorted.bg"
        echo "Sorting the cleaned bedgraph..." >> $logfile 2>&1
        sortBed -i $bedgraph_output_clean > $bedgraph_sort_output
        """


# skipping step 7


rule L1_loci_filter:
    input:
        step6=rules.read_filter.output,
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_consistent_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_threshold_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.sorted.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/raw_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/filtered_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/active_coverage_for_weighted_avg.bed",
    conda:
        "line_expression_lrs.yaml"
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/L1_loci_filter.log",
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        L1_ref_type="active",
    shell:
        """
        logfile="../{params.sample}/log/L1_loci_filter.log"

        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} {params.L1_ref_type} > $logfile 2>&1
        """


rule final_map_qc_LRS:
    input:
        step8=rules.L1_loci_filter.output,
        bam_input=rules.read_filter.output[5],
    output:
        directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}"
        ),
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/09_map_qc_LRS.log",
    shell:  # TODO add threads
        """
        logfile="../../../log/09_map_qc_LRS.log"

        cd results/LINE-Expression-LRS/{params.sample}/d_LINE_quantification/active/read_filter/
        longreadsum bam -i $(basename {input.bam_input}) -o {params.sample} > $logfile 2>&1
        """


rule normalization_wgt_avg:
    input:
        step9=rules.final_map_qc_LRS.output,
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/normalized_active_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/coverage_weighted_avg.bed",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        L1_ref_type="active",
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/normalization_wgt_avg.log",
    shell:
        """
        logfile="../{params.sample}/log/normalization_wgt_avg.log"

        cd results/LINE-Expression-LRS/scripts
        ./$(basename {input.script}) {params.sample} {params.L1_ref_type} > $logfile 2>&1

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
                rules.preprocess_mapping.output,
                zip,
                sample=ss["sample"],
                libtype=ss.libtype,
            )


rule lrs_done:
    input:
        get_lrs_output,
