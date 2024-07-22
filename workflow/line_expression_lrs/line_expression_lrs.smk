
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
        fasta="results/LINE-Expression-LRS/{sample}_{libtype}/reads.fa",
        fasta1kb="results/LINE-Expression-LRS/{sample}_{libtype}/reads_cDNA_1kb.fa",
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


rule get_l1hs_hmm:
    input:
        download_complete=rules.download_line_expresssion_lrs.output[0],
    output:
        hmm="resources/L1HS_5end.hmm",
    log:
        "resources/L1HS_5end_dfam_query.log",
    conda:
        "line_expression_lrs.yaml"
    shell:
        """
        curl -s https://dfam.org/api/families/DF000000226/hmm?format=hmm > {output.hmm} 2> {log}
        """


rule L1_detection:
    input:
        fa=rules.preprocess_input.output.fasta1kb,
        lib=rules.get_l1hs_hmm.output.hmm,
    output:
        rmsk=multiext(
            rules.preprocess_input.output.fasta1kb,
            ".align",
            ".masked",
            ".ori.out",
            ".out",
            ".out.xm",
            ".tbl",
        ),
        ids="results/LINE-Expression-LRS/{sample}_{libtype}/l1hs_reads.txt",
        fa="results/LINE-Expression-LRS/{sample}_{libtype}/l1hs_reads.fa",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    threads: 32
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/L1_detection.log",
    shell:
        """
        exec &>> {log}

        echo "Running RepeatMasker at $(date)..."
        ## FLAGS
        # -div = max divergence
        # -lib = library of sequences to search for
        # -no_is = skips bacterial insertion element check
        # -nolow = does not mask low complexity DNA or simple repeats
        # -norna = no RNA repeats. When interested in small RNA genes, you should use the -norna option that leaves these sequences unmasked, while still masking SINEs.
        # -a = shows the alignments in a .align output file
        # -u = creates an untouched annotation file besides the manipulated file
        # -xsmall = returns repetitive regions in lowercase (rest capitals) rather than masked
        # -xm = creates an additional output file in cross_match format (for parsing)
        # -e hmmer = use HMMER for engine
        # -s = slow search / -qq = Rush job; about 10% less sensitive,

        RepeatMasker -pa {threads} -lib {input.lib} -no_is -nolow -norna -a -u -xsmall -xm -div 10 -e hmmer -qq {input.fa}

        echo "Finished RepeatMasker at $(date)!"

        # Post-RepeatMasker Filtering by 10% Divergence
        echo "Filtering RepeatMasker output by 10% divergence at $(date)..."

        awk '$0 ~ /L1HS/ && $2 <= 10 {{print $5}}' {output.rmsk[3]} > {output.ids}

        echo "Generating a new FASTA file of reads with less than 10% diverged LINE/L1 elements at $(date)..."

        seqtk subseq {input.fa} {output.ids} > {output.fa}

        echo "Completed processing RepeatMasker output and generated a new FASTA file of reads with less than 10% diverged LINE/L1 elements!"
        """


rule map_hg38:
    input:
        fa=rules.L1_detection.output.fa,
        ref=remote_or_local(config["genome_fa"]),
        gen=remote_or_local(config["gencode_gtf"]),
    output:
        bam="results/LINE-Expression-LRS/{sample}_{libtype}/hg38_mapped.sorted.bam",
        bai="results/LINE-Expression-LRS/{sample}_{libtype}/hg38_mapped.sorted.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
    threads: 8
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/map_hg38.log",
    shell:
        """
        exec &>> {log}

        echo "Mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome..."

        minimap2 -ax splice --junc-bed {input.gen} -uf --secondary=no -k14 -t {threads} {input.ref} {input.fa} | \
            samtools view -h -@ {threads} - | \
            samtools sort -@ {threads} - > {output.bam}

        samtools index {output.bam}

        echo "Finished mapping reads with < 10% LINE/L1 elements to the hg38 Reference Genome"
        """


rule map_qc_LRS:
    input:
        bam_input=rules.map_hg38.output.bam,
    output:
        directory("results/LINE-Expression-LRS/{sample}_{libtype}/hg38_mapping_LRS"),
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
        L1_ref_input="resources/LINE-Expression-LRS/references/L1Base2_filtered/active_filtered.bed",
        bam=rules.map_hg38.output.bam,
        map_qc_LRS=rules.map_qc_LRS.output[0],
    output:
        bedgraph_output_clean="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_bedgraph_clean.bg",
        bedgraph_sort_output="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_bedgraph_sorted.bg",
        bedgraph_output="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_bedgraph.bg",
        L1_regions_reads="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_L1_regions_reads.bam",
        read_filter_bam="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_read_filter_passed.bam",
        sorted_read_filter_bam="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam",
        index="results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/{sample}_{libtype}_read_filter_passed.sorted_position.bam.bai",
    conda:
        "line_expression_lrs.yaml"
    params:
        sample=lambda wc: wc.sample + "_" + wc.libtype,
        L1_ref_type="active",
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/read_filter.log",
    shell:
        """
        exec &>> {log}

        echo "Read Filter on the {params.L1_ref_type} Reference L1 Regions"

        # Output Files
        L1_regions_reads={output.L1_regions_reads}
        echo "Generating filtered BAM file with reads only located within the L1 reference regions..."
        samtools view -b -L {input.L1_ref_input} -o "$L1_regions_reads" {input.bam}

        read_filter_bam={output.read_filter_bam}
        echo "Removing reads with less than 90% of the read maps to the L1 reference regions"
        bedtools intersect -a "$L1_regions_reads" -b "$L1_regions_reads" -f 0.9 > "$read_filter_bam"

        sorted_read_filter_bam={output.sorted_read_filter_bam}
        echo "Sorting and Indexing the resulting Read Filter BAM file..."
        samtools sort "$read_filter_bam" -o "$sorted_read_filter_bam"
        samtools index "$sorted_read_filter_bam"

        # Generate the bedgraph for the new BAM file that passed the Read Filter
        bedgraph_output={output.bedgraph_output}
        echo "Generating the bedgraph..."
        bedtools genomecov -ibam "$sorted_read_filter_bam" -bga -split > "$bedgraph_output"

        bedgraph_output_clean={output.bedgraph_output_clean}
        echo "Cleaning the bedgraph..."
        grep -v 'fix\|alt\|random\|[(]\|Un' $bedgraph_output > $bedgraph_output_clean

        bedgraph_sort_output={output.bedgraph_sort_output}
        echo "Sorting the cleaned bedgraph..."
        sortBed -i $bedgraph_output_clean > $bedgraph_sort_output
        """


rule final_map_qc_LRS:
    input:
        bam_input=rules.read_filter.output.sorted_read_filter_bam,
    output:
        directory(
            "results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/read_filter_LRS"
        ),
        "results/LINE-Expression-LRS/{sample}_{libtype}/read_filter/read_filter_LRS/bam_summary.txt",
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
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/{sample}_{libtype}_consistent_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/{sample}_{libtype}_threshold_passed_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/{sample}_{libtype}_regions_for_coverage.sorted.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/raw_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/filtered_coverage_values_mean.txt",
        "results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter/active_coverage_for_weighted_avg.bed",
        directory("results/LINE-Expression-LRS/{sample}_{libtype}/L1_loci_filter"),
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
        "results/LINE-Expression-LRS/{sample}_{libtype}/normalized_active_regions.bed",
        "results/LINE-Expression-LRS/{sample}_{libtype}/coverage_weighted_avg.bed",
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


rule LRS_reads_report_preprocess:
    input:
        original_fq=get_lrs_fq,
        preprocess_input=rules.preprocess_input.output.fasta1kb,
        L1_detection=rules.L1_detection.output.fa,
        map_hg38=rules.map_hg38.output.bam,
        read_filter=rules.read_filter.output.sorted_read_filter_bam,
        L1_loci_filter=rules.normalization_wgt_avg.input.input_ref_cov,
        normalization_wgt_avg=rules.normalization_wgt_avg.output[0],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/LRS_reads_report_preprocess.csv",
    conda:
        "line_expression_lrs.yaml"
    log:
        "results/LINE-Expression-LRS/{sample}_{libtype}/log/LRS_reads_report_preprocess.log",
    shell:
        """
        exec &>> {log}
        # get read counts for each step's inputs/output

        if [[ {input.original_fq} == *.gz ]]; then
            original_fq_read_count=$(zgrep -c '^@' {input.original_fq})
        else
            original_fq_read_count=$(grep -c '^@' {input.original_fq})
        fi
        preprocess_input_read_count=$(grep -c '^>' {input.preprocess_input})
        L1_detection_read_count=$(grep -c '^>' {input.L1_detection})
        map_hg38_read_count=$(samtools view -c {input.map_hg38})
        read_filter_read_count=$(samtools view -c {input.read_filter})
        L1_loci_filter_read_count=$(wc -l {input.L1_loci_filter} | awk '{{print $1}}')
        normalization_wgt_avg_read_count=$(wc -l {input.normalization_wgt_avg} | awk '{{print $1}}')

        echo "original_fq,preprocess_input,L1_detection,map_hg38,read_filter,L1_loci_filter,normalization_wgt_avg" > {output[0]}
        echo "$original_fq_read_count,$preprocess_input_read_count,$L1_detection_read_count,$map_hg38_read_count,$read_filter_read_count,$L1_loci_filter_read_count,$normalization_wgt_avg_read_count" >> {output[0]}

        """


rule LRS_reads_report:
    input:
        reads=rules.LRS_reads_report_preprocess.output[0],
    output:
        "results/LINE-Expression-LRS/{sample}_{libtype}/LRS_reads_report.ipynb",
    conda:
        "line_expression_lrs.yaml"
    log:
        notebook="results/LINE-Expression-LRS/{sample}_{libtype}/LRS_reads_report.ipynb",
    notebook:
        "LRS_reads_report.py.ipynb"


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
                rules.LRS_reads_report.output,
                zip,
                sample=ss["sample"],
                libtype=ss.libtype,
            )


rule lrs_done:
    input:
        get_lrs_output,
