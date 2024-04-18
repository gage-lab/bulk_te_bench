#!/bin/bash

bam_files=(
    "../illumina_files/SGNex_MCF7_Illumina_replicate2_run1/SGNex_MCF7_Illumina_replicate2_run1.bam"
    "../illumina_files/SGNex_MCF7_Illumina_replicate3_run1/SGNex_MCF7_Illumina_replicate3_run1.bam"
    "../illumina_files/SGNex_MCF7_Illumina_replicate4_run1/SGNex_MCF7_Illumina_replicate4_run1.bam"
)


for i in "${bam_files[@]}"; do
	echo "Processing file: $i"
    output_file="${i%.bam}_unique_reads.txt"
    samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c > "$output_file"
    python 02_bedtools_gtf.py "$i"
done
