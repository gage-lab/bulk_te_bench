#!/bin/bash


## SGNex Files

bam_files=(
    "../illumina_files/SGNex_MCF7_Illumina_replicate2_run1/SGNex_MCF7_Illumina_replicate2_run1.bam"
    "../illumina_files/SGNex_MCF7_Illumina_replicate3_run1/SGNex_MCF7_Illumina_replicate3_run1.bam"
    "../illumina_files/SGNex_MCF7_Illumina_replicate4_run1/SGNex_MCF7_Illumina_replicate4_run1.bam"
)

for i in "${bam_files[@]}"; do
	echo "Processing file: $i"
    output_file="${i%.bam}_unique_reads.txt"
    #samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c | awk '$1 == 2' | awk '{print $2}' > "$output_file"
    #python 02_bedtools_gtf.py "$i"
done


## Mike's Files

output_dir="../illumina_files/mikes_files"

### Genome aligned

bam_files=(
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate2_run1/Aligned.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate3_run1/Aligned.out.bam"
    "../illumina_files/mikes_files/star_align/SGNex_MCF7_Illumina_replicate4_run1/Aligned.out.bam"
)


for i in "${bam_files[@]}"; do
    echo "Processing file: $i"
    replicate_run=$(basename $(dirname "$i"))
    prefix="${output_dir}/genome/${replicate_run}"
    output_file="${prefix}/unique_reads.txt"
    #samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c | awk '$1 == 2'| awk '{print $2}' > "$output_file"
    python 02_bedtools_gtf.py "$i" --outdir "${prefix}" --chr
done

### Transcriptome aligned

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
    #samtools sort -n "$i" | samtools view -f 3 | awk '{print $1}' | uniq -c | awk '$1 == 2' | awk '{print $2}' > "$output_file"
done
