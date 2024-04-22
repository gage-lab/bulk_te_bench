#!/bin/bash

# bam_files=(
#     "../illumina_files/SGNex_MCF7_Illumina_replicate2_run1/SGNex_MCF7_Illumina_replicate2_run1.bam"
#     "../illumina_files/SGNex_MCF7_Illumina_replicate3_run1/SGNex_MCF7_Illumina_replicate3_run1.bam"
#     "../illumina_files/SGNex_MCF7_Illumina_replicate4_run1/SGNex_MCF7_Illumina_replicate4_run1.bam"
# )
#
#
# for i in "${bam_files[@]}"; do
# 	echo "Processing file: $i"
#     output_file="${i%.bam}_unique_reads.txt"
#     samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c > "$output_file"
#     python 02_bedtools_gtf.py "$i"
# done


bam_files=(
    #"../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate2_run1/Aligned.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate3_run1/Aligned.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate4_run1/Aligned.out.bam"
)

output_dir="../illumina_files/mikes_files"

for i in "${bam_files[@]}"; do
    echo "Processing file: $i"
    replicate_run=$(basename $(dirname "$i"))
    prefix="${output_dir}/genome/${replicate_run}"
    output_file="${prefix}/unique_reads.txt"
    #samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c > "$output_file"
    python 02_bedtools_gtf.py "$i" --outdir "${prefix}"
done


bam_files=(
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate2_run1/Aligned.toTranscriptome.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate3_run1/Aligned.toTranscriptome.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate4_run1/Aligned.toTranscriptome.out.bam"
)

for i in "${bam_files[@]}"; do
    echo "Processing file: $i"
    replicate_run=$(basename $(dirname "$i"))
    prefix="${output_dir}/txome/${replicate_run}"
    output_file="${prefix}/unique_reads.txt"
    #samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c > "$output_file"
    #python 02_bedtools_gtf.py "$i" "--outdir ${prefix}"
done
