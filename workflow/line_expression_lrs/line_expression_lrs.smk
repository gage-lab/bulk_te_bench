
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
        "resources/LINE-Expression-LRS/references/L1Base2_filtered/active_filtered.bed",
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
        mkdir -p results/LINE-Expression-LRS/{params.sample}/a_dataset
        mkdir -p results/LINE-Expression-LRS/{params.sample}/b_repeat_masker_process
        mkdir -p results/LINE-Expression-LRS/{params.sample}/c_hg38_mapping_LRS
        mkdir -p results/LINE-Expression-LRS/{params.sample}/d_LINE_quantification

        echo "All folders have been made!" >> {log}


        # Input File Preparation

        if [ {params.filetype} == "FASTQ" ]; then
            echo "Input file is:   FASTQ" >> {log}
            echo "Converting to FASTA..." >> {log}
            seqtk seq -a {input.fastq} > {output.fasta} #check good


            echo "Filtering out reads less than 1kb..." >> {log}
            seqtk seq -L 1000 {output.fasta} > {output.fasta1kb} # #TODO: MAKE THIS TEMP OUTPUT

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
        samtools fasta {output.sam} -F 2308 -@ {threads} > {output.fa} 2>{log}

        echo "Mapping to Custom LINE Reference Library complete!" >> {log}

        """


rule L1_detection:
    input:
        input_file=rules.preprocess_mapping.output.fa,
        fasta1kb=rules.preprocess_input.output.fasta1kb,
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

        echo "Running RepeatMasker..." >> {log} 2>&1
        RepeatMasker -pa 6 -dir "results/LINE-Expression-LRS/{params.sample}/b_repeat_masker_process" -nolow -norna -div 10 -species human -no_is -a -u -xsmall -xm {input.input_file} >> {log} 2>&1

        echo "Finished RepeatMasker!" >> {log} 2>&1

        # Post-RepeatMasker Filtering by 10% Divergence

        RM_input_file={output[3]}
        output_file={output[7]}

        readIDs_lines=()

        while read -r line; do
            if [[ $line == *"LINE/L1"* ]]; then
                fields=($line)

                if (( $(echo "${{fields[1]}} <= 10" | bc -l) )); then
                    echo "$line" >> "$output_file"
                    readIDs_lines+=("${{fields[4]}}")
                fi
            fi
        done < "$RM_input_file"



        # Get ReadIDs and Generate the new FASTA file of these ReadIDs
        read_id_file={output[8]}

        echo "Gathering Final ReadIDs of less than 10% diverged LINE/L1 elements..." >> {log} 2>&1

        for item in "${{readIDs_lines[@]}}"; do
            echo "$item" >> "$read_id_file"
        done


        seqtk subseq {input.fasta1kb} $read_id_file > {output[6]}

        echo "Completed processing RepeatMasker output and generated a new FASTA file of reads with less than 10% diverged LINE/L1 elements!" >> {log} 2>&1

        """


rule map_hg38:
    input:
        fasta_input=rules.L1_detection.output[6],
        ref=remote_or_local(config["genome_fa"]),
        gen=remote_or_local(config["gencode_gtf"]),
    output:
        sam="results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sam",
        bam="results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.bam",
        sorted_bam="results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam",
        index="results/LINE-Expression-LRS/{sample}_{libtype}/c_hg38_mapping_LRS/{sample}_{libtype}_hg38_mapped.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    threads: 8
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/map_hg38.log",
    shell:
        """
        echo "Mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome..." > {log} 2>&1

        minimap2 -ax splice --junc-bed {input.gen} -uf --secondary=no -k14 -t {threads} {input.ref} {input.fasta_input} -o {output.sam}  >> {log} 2>&1
        samtools view -Sb -o {output.bam} {output.sam}  >> {log} 2>&1
        samtools sort {output.bam} -@ {threads} -o {output.sorted_bam}  >> {log} 2>&1
        samtools index {output.sorted_bam}  >> {log} 2>&1

        echo "Finished mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome" >> {log} 2>&1
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
        longreadsum bam -i {input.bam_input} -o {output[0]} > {log} 2>&1

        """


### after this we need a new parameter: active, inactive, or ORF2. Moving forwards with active for now for now
rule read_filter:  # TODO turn every instance of "active" into a wc
    input:
        step5=rules.map_qc_LRS.output,
        L1_ref_input="resources/LINE-Expression-LRS/references/L1Base2_filtered/active_filtered.bed",
        sample_bam_input=rules.map_hg38.output.sorted_bam,
    output:
        bedgraph_output_clean="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph_clean.bg",
        bedgraph_sort_output="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph_sorted.bg",
        bedgraph_output="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_bedgraph.bg",
        L1_regions_reads="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_L1_regions_reads.bam",
        read_filter_bam="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.bam",
        sorted_read_filter_bam="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam",
        index="results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        L1_ref_type="active",
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/read_filter.log",
    shell:
        """
        echo "Read Filter on the {params.L1_ref_type} Reference L1 Regions" > {log} 2>&1

        mkdir -p results/LINE-Expression-LRS/{params.sample}/d_LINE_quantification/{params.L1_ref_type}
        mkdir -p results/LINE-Expression-LRS/{params.sample}/d_LINE_quantification/{params.L1_ref_type}/read_filter

        # Output Files
        L1_regions_reads={output.L1_regions_reads}
        echo "Generating filtered BAM file with reads only located within the L1 reference regions..." >> {log} 2>&1
        samtools view -b -L {input.L1_ref_input} -o "$L1_regions_reads" {input.sample_bam_input}

        read_filter_bam={output.read_filter_bam}
        echo "Removing reads with less than 90% of the read maps to the L1 reference regions" >> {log} 2>&1
        bedtools intersect -a "$L1_regions_reads" -b "$L1_regions_reads" -f 0.9 > "$read_filter_bam"

        sorted_read_filter_bam={output.sorted_read_filter_bam}
        echo "Sorting and Indexing the resulting Read Filter BAM file..." >> {log} 2>&1
        samtools sort "$read_filter_bam" -o "$sorted_read_filter_bam"
        samtools index "$sorted_read_filter_bam"

        # Generate the bedgraph for the new BAM file that passed the Read Filter
        bedgraph_output={output.bedgraph_output}
        echo "Generating the bedgraph..." >> {log} 2>&1
        bedtools genomecov -ibam "$sorted_read_filter_bam" -bga -split > "$bedgraph_output"

        bedgraph_output_clean={output.bedgraph_output_clean}
        echo "Cleaning the bedgraph..." >> {log} 2>&1
        grep -v 'fix\|alt\|random\|[(]\|Un' $bedgraph_output > $bedgraph_output_clean

        bedgraph_sort_output={output.bedgraph_sort_output}
        echo "Sorting the cleaned bedgraph..." >> {log} 2>&1
        sortBed -i $bedgraph_output_clean > $bedgraph_sort_output
        """


# skipping step 7


rule final_map_qc_LRS:
    input:
        bam_input=rules.read_filter.output.sorted_read_filter_bam,
    output:
        directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}"
        ),
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/read_filter/{sample}_{libtype}/bam_summary.txt",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/final_map_qc_LRS.log",
    shell:
        """
        longreadsum bam -i {input.bam_input} -o {output[0]} > {log} 2>&1
        """


rule L1_loci_filter:
    input:
        sorted_output_bam=rules.read_filter.output.sorted_read_filter_bam,
        bedgraph_sort_output=rules.read_filter.output.bedgraph_sort_output,
        L1_ref_regions="resources/LINE-Expression-LRS/references/L1Base2_filtered/active_filtered.bed",
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_consistent_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_threshold_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.sorted.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/raw_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/filtered_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter/active_coverage_for_weighted_avg.bed",
        directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/d_LINE_quantification/active/L1_loci_filter"
        ),
    conda:
        "line_expression_lrs.yaml"
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/L1_loci_filter.log",
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        L1_ref_type="active",
    script:
        "L1_loci_filter.sh"


rule normalization_wgt_avg:
    input:
        summary=rules.final_map_qc_LRS.output[1],
        input_ref_cov=rules.L1_loci_filter.output[6],
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
        echo "Normalization by Total Number of Reads..." >> {log}
        qc_report_input={input.summary}

        # function to calculate coverage

        # Extract the number of read value from summary.txt
        value=$(awk -F '\t' 'NR==1 {{print $2}}' "$qc_report_input")

        #value=$(awk '{{print $5}}' $qc_report_input)
        echo "Number of reads: $value" >> {log}

        input_ref_cov={input.input_ref_cov}


        # output files
        normalized_ref_regions={output[0]}

        awk -v divisor="$value" -v OFS="\t" '{{$NF = (divisor != 0) ? $NF / divisor : 0; print}}' "$input_ref_cov" > "$normalized_ref_regions"

        echo "Normalization Complete!" >> {log}


        # Part 2: Weighted Average Calculation
        echo "Weighted Average Calculation..." >> {log}

        # output file
        weighted_average_cov={output[1]}


        weighted_sum=0
        total_weight=0

        while IFS=$'\t' read -r line || [[ -n "$line" ]]; do
            elements=($line)
            start=${{elements[1]}}
            end=${{elements[2]}}
            coverage=$(echo "${{elements[-1]}}" | awk '{{print $NF}}')
            region_size=$((end - start + 1))
            weighted_sum=$(awk "BEGIN {{print $weighted_sum + ($coverage * $region_size)}}")
            total_weight=$((total_weight + region_size))
        done < "$input_ref_cov"

        if [ "$total_weight" -ne 0 ]; then
            weighted_average=$(awk "BEGIN {{print $weighted_sum / $total_weight}}")
        else
            weighted_average=0
        fi

        # Write the Sample Name and Coverage Values to the new file
        echo -e "Sample Name\tWeighted Average" > "$weighted_average_cov"
        echo -e "{params.sample}\t$weighted_average" >> "$weighted_average_cov"

        echo "Calculations for {params.sample}, over the {params.L1_ref_type} regions is complete!" >> {log}
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
                rules.L1_loci_filter.output,
                zip,
                sample=ss["sample"],
                libtype=ss.libtype,
            )


rule lrs_done:
    input:
        get_lrs_output,
