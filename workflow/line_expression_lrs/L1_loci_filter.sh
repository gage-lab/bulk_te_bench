#!/usr/bin/env bash

echo "L1 Loci Filter on the ${params[L1_ref_type]} Reference L1 Regions" >> "${snakemake_snakemake_log}"

sorted_output_bam="${snakemake_input[sorted_output_bam]}"
L1_ref_regions="${snakemake_input[L1_ref_regions]}"
bedgraph_sort_output="${snakemake_input[bedgraph_sort_output]}"

mkdir -p "results/LINE-Expression-LRS/a_DNA/d_LINE_quantification/active/L1_loci_filter"


###########
# Check 1 #
###########
# Checking the read start or end positions (taking into account strandness) and filtering regions with inconsistent read start positions
echo "Checking the read start or end positions (taking into account strandness) and filtering regions with inconsistent read start positions" >> "${snakemake_log}"

# Set the output file paths
CONSISTENT_REGIONS_FILE="${snakemake_output[0]}"

# Create the output files (or clear their contents if they already exist)
> "$CONSISTENT_REGIONS_FILE"

# Read from the original reference regions file ($OG_REF_REGIONS) instead
while IFS=$'\t' read -r chrom start end name score strand; do
    TEMP_FILE=$(mktemp)
    samtools view -b "$sorted_output_bam" "$chrom:$start-$end" | bedtools bamtobed -i - > "$TEMP_FILE" 2> "${snakemake_log}"

    # Check if the temporary file is empty
    if [ ! -s "$TEMP_FILE" ]; then
        continue
    fi

    # Determine whether to check starting position or end position based on strand
    if [[ "$strand" == "+" ]]; then
        position_col=2
        position_label="start"
    else
        position_col=3
        position_label="end"
    fi

    # Extract the end positions of the reads within the region
    end_positions=$(awk -v pos_col="$position_col" '{print $pos_col}' "$TEMP_FILE")

    # Check if end positions of reads are consistent within 100 bps of each other
    consistent_count=$(echo "$end_positions" | awk -v diff_limit=100 '{
        prev_pos=0;
        count=0;
        for (i=1; i<=NF; i++) {
            if (prev_pos != 0 && ($i - prev_pos) > diff_limit) {
                exit 1;
            }
            prev_pos = $i;
            count++;
        }
    }') >> "${snakemake_log}"

    if [ -z "$consistent_count" ]; then
        # Output the region to the consistent regions file
        echo -e "$chrom\t$start\t$end\t$name\t$score\t$strand" >> "$CONSISTENT_REGIONS_FILE"
    fi

    # Remove the temporary file
    rm "$TEMP_FILE"

done < "$L1_ref_regions"

echo "Finished Check #1" >> "${snakemake_log}"


###########
# Check 2 #
###########
echo "Checking if the starting position falls within the 1.5kb window between the average consistent starting position and the reference starting position (taking into account strandness)." >> "${snakemake_log}"


# Set the output file paths
OUTPUT_FILE_threshold="${snakemake_output[2]}"
UNDER1500_FILE="${snakemake_output[1]}"


# Create the output files (or clear their contents if they already exist)
> "$OUTPUT_FILE_threshold"
> "$UNDER1500_FILE"

# Process positive strand regions
awk '$6 == "+" {print}' "$L1_ref_regions" | while IFS=$'\t' read -r chrom start end name score strand; do
    TEMP_FILE=$(mktemp)
    samtools view -b "$sorted_output_bam" "$chrom:$start-$end" | bedtools bamtobed -i - > "$TEMP_FILE" 2> "${snakemake_log}"

    # Check if the temporary file is empty
    if [ ! -s "$TEMP_FILE" ]; then
        continue
    fi

    position_col=2
    position_label="start"

    # Extract the positions of the reads within the region
    positions=$(cut -f "$position_col" "$TEMP_FILE")

    # Calculate the differences
    differences=()
    for pos in $positions; do
        diff=$((pos - start))
        differences+=("$diff")
    done

    # Find the mode of the differences
    mode_diff=$(printf '%s\n' "${differences[@]}" | awk '{a[$1]++}END{for(i in a){if(a[i]>max){max=a[i];n=i}}}END{print n}')

    # Check if the mode is negative and make it positive
    if [ $mode_diff -lt 0 ]; then
        mode_diff=$((-$mode_diff))
    fi

    # Check if the mode difference is below 1500
    if [ $mode_diff -lt 1500 ]; then
        echo "Region: ${chrom}_${start}_${end} (Strand: $strand)" >> "${snakemake_log}"
        echo "Mode Difference: $mode_diff" >> "${snakemake_log}"

        # Append the region to the under 1500 file
        echo -e "${chrom}\t${start}\t${end}\t${name}\t${score}\t${strand}" >> "$UNDER1500_FILE"
    fi

    # Append the extracted reads to the output file
    cat "$TEMP_FILE" >> "$OUTPUT_FILE_threshold"

    # Remove the temporary file
    rm "$TEMP_FILE"

done

# Process negative strand regions
# the input bam file for this part is under the variable:  $SORTED_OUTPUT_BAM
# the input reference regions are under $OG_REF_REGIONS


awk '$6 == "-" {print}' "$L1_ref_regions" | while IFS=$'\t' read -r chrom start end name score strand; do
    TEMP_FILE=$(mktemp)
    samtools view -b "$sorted_output_bam" "$chrom:$start-$end" | bedtools bamtobed -i - > "$TEMP_FILE" 2> "${snakemake_log}"

    # Check if the temporary file is empty
    if [ ! -s "$TEMP_FILE" ]; then
        continue
    fi

    position_col=3
    position_label="end"

    # Extract the positions of the reads within the region
    positions=$(cut -f "$position_col" "$TEMP_FILE")

    # Perform your analysis on the read positions here
    # Replace the following echo statements with your desired snakemake_logic

    # Calculate the differences
    differences=()
    for pos in $positions; do
        diff=$((end - pos))
        differences+=("$diff")
    done

    # Find the mode of the differences
    mode_diff=$(printf '%s\n' "${differences[@]}" | awk '{a[$1]++}END{for(i in a){if(a[i]>max){max=a[i];n=i}}}END{print n}')

    # Check if the mode is negative and make it positive
    if [ $mode_diff -lt 0 ]; then
        mode_diff=$((-$mode_diff))
    fi

    # Check if the mode difference is below 1500
    if [ $mode_diff -lt 1500 ]; then
        echo "Region: ${chrom}_${start}_${end} (Strand: $strand)" >> "${snakemake_log}"
        echo "Mode Difference: $mode_diff" >> "${snakemake_log}"

        # Append the region to the under 1500 file
        echo -e "${chrom}\t${start}\t${end}\t${name}\t${score}\t${strand}" >> "$UNDER1500_FILE"
    fi

    # Append the extracted reads to the output file
    cat "$TEMP_FILE" >> "$OUTPUT_FILE_threshold"

    # Remove the temporary file
    rm "$TEMP_FILE"

done



# sort the regions to look over
coverage_regions="${snakemake_output[3]}"
bedtools sort -i $UNDER1500_FILE > $coverage_regions

echo "Finished Check #2" >> "${snakemake_log}"


###########
# Check 3 #
###########
echo "Calculate coverage over the remaining L1 regions and filter those with less than 2 reads" >> "${snakemake_log}"


coverage_output_mean="${snakemake_output[4]}"
bedtools map -a $coverage_regions -b $bedgraph_sort_output -c 4 -o mean -null 0 > $coverage_output_mean 2> "${snakemake_log}"
echo "calculated coverage; by MEAN" >> "${snakemake_log}"

# Filter regions with a value less than 3 in the last column
filtered_coverage_output_mean="${snakemake_output[5]}"
awk '$NF >= 3' $coverage_output_mean > $filtered_coverage_output_mean
echo "Filtered regions with less than 2 reads. " >> "${snakemake_log}"


# Replace regions with less than 2 reads with 0
FINAL_COV_OUTPUT="${snakemake_output[6]}"

if [[ -s $filtered_coverage_output_mean ]]; then
    awk 'NR==FNR{regions[$1,$2,$3]=$0; next} {if (($1,$2,$3) in regions) print regions[$1,$2,$3]; else print $0}' $filtered_coverage_output_mean $OG_REF_REGIONS > $FINAL_COV_OUTPUT
else
    cp $L1_ref_regions $FINAL_COV_OUTPUT
fi

echo "Calculated the "${params[L1_ref_type]}" regions coverage values" >> "${snakemake_log}"
