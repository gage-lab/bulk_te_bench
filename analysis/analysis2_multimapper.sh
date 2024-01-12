#!/bin/bash

# run code in "analysis" directory

# Directories containing BAM files
directories=("longread_files" "illumina_files")

# Iterate over each directory
for input_directory in "${directories[@]}"; do
    # Iterate over BAM files in the directory
    for bam_file in "$input_directory"/*/*.bam; do
        if [ -e "$bam_file" ]; then
            echo $bam_file

            # Extract the filename without extension
            bam_id=$(basename "$bam_file" .bam)

            # Step 1: Convert BAM to TSV
            samtools view $bam_file | awk -F'\t' '{score=""; for(i=3; i<=NF; i++) { if ($i ~ /^AS:i:/) { score=substr($i, 6); break } } print $1 "\t" score}' >  "$input_directory/$bam_id/$bam_id.tsv"

            # Step 2: Sort TSV by Read ID and Score
            sort -k1,1 -k2,2nr -o "$input_directory/$bam_id/$bam_id.sorted.tsv" "$input_directory/$bam_id/$bam_id.tsv"

            # Step 3: Count Alignments with Top Score and Create Lookup Table
            awk -F'\t' -v OFS='\t' '{
                if ($2 == prev_score) {
                    count++;
                } else {
                    if (NR > 1) {
                        print prev_readID, count;
                    }
                    prev_readID = $1;
                    prev_score = $2;
                    count = 1;
                }
            }
            END {
                print prev_readID, count;
            }' "$input_directory/$bam_id/$bam_id.sorted.tsv" > "$input_directory/$bam_id/${bam_id}_lookup_table.txt"
        fi

    done
done
