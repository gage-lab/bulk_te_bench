This analysis folder contains several analyses for looking at L1HS expression in the [Singapore Nanopore Expression Project](https://github.com/GoekeLab/sg-nex-data).

[Slides](https://docs.google.com/presentation/d/1dDzzKKmCZ8zZJ1IdnnbIdjcjSQFwULDBAFj0lTNhwRA/edit#slide=id.g28e8f426d9c_0_0)
# Structure
Data used for analyses is generated in `preview.ipnyb` and stored in respective directories, or the `data/` directory.

## preview.ipynb
Contains test code and analysis code. When all analyses are good and done, will remove those test sections from code.

Analysis sections:
1.  "Split up L1HS into truncated, full length, intronic, intergenic"
2.  "Calculate total signal per bigwig"

# Getting data and generating heatmaps
Purpose: Download relevant [SG-NEx data](https://github.com/GoekeLab/sg-nex-data) and mae
## Files
1. `main.py`
- Download long_read bigwigs from AWS and generates heatmaps
2. `main_illumina.py`
- Download illumina bams from AWS, turns into bigwigs, and generates heatmaps
3. `preview.ipynb`: "Split up L1HS into truncated, full length, intronic, intergenic"
- Splits up L1HS bed files for heatmaps. Store output in `data/`
4. `generate_heatmap.sh`
- Generates heatmaps using bed files and inputted bigwig file. Outputs in `heatmaps/`
5. `download_long.py`
- Downloads long read bam files for future analyses

# Signals per bigwig file

After getting bigwigs, create plots demonstrating signal between different groups. This is kind of preliminary.

## Output
1. Line plots with signal values for each bigwig
2. Strip plots of signals intersecting with bed file regions

## Files and relevant functions
`preview.ipynb`: "Calculate total signal per bigwig"
Functions -- for all bigwig files in a directory, caclulate the the sum of all values in the file that intersect with a bed file
- `calculate_total_signal(bigwig_path, bed_path)`
- `calculate_average_signal(bigwig_path, bed_path)`
** bw.stats() returns mean unless otherwise states

# Lookup table to find multimapper sites for all bw files
How many multimappers occur at L1 sites in illumina vs nanopore data? Relative to the rest of the genome?

## Files

`analysis2_multimapper.sh`:
1. bash code that iterates over both directories (long read and illumina) and gets all bam files
2. Turns bam files into TSVs of read_id, alignment score
3. sort TSV by read_ID then by score
4. count number of alignments with top score
5. save to lookup table of read_id, counts

## Next steps

1. Pick reads of interest
- Use bedtools to query bam file to get all read_ids based off of bed file with loci of interests
- Use that list of read_ids to filter the look_up table per bam file
2. Use tables to compare counts across genome for reads of interest
- check all the reads that align to l1hs seqs, then randomly sample the same # of reads 100(0) times (background distribution), and then test if statistically higher multi-mapping for that location
- take mean for each distribution --> jitter plot
- can start with MCF-7 direct-cDNA, can compare to short reads too



# Aggregate data for replicates across runs

## Files


## new idea
Look in pape/github how they split up replicates + runs (https://github.com/GoekeLab/sg-nex-data/issues/55)
1. merge bam files
2. re-sort
3. make bigwig (copy process for illumina bws)


This post warns against this -- but weirdly tough to do.

Why is it tough?
- No smooth tool for doing so (UCSC's `bigWigMerge` returns bedgraph)
- The bws have lots of alternatives in chromosomes so that must be settled
- For some reason there is an extra col in final bedgraph (last two cols are filled with [[0, 1]])

## Code looks something like this

```
bigWigToBedGraph RepX_RunA.bw run1.bedGraph
bigWigToBedGraph RepX_RunB.bw run2.bedGraph

grep -wFf chrom.sizes run1.bedGraph > run1_fil.bedGraph
grep -wFf chrom.sizes run2.bedGraph > run2_fil.bedGraph


bedtools unionbedg -i run1_fil.bedGraph run2_fil.bedGraph > merged.bedGraph
awk '{print $1, $2, $3, $4}' merged.bedGraph > merged_fixed.bedGraph
bedGraphToBigWig merged_fixed.bedGraph hg38.chrom.sizes merged.bigwig
bash generate_heatmap.sh merged.bigwig

```

# WIP - How often do we observe reads at expected loci? Full-length vs non full-length intergenic L1HS?

Utilize `calculate_total_signal(bigwig_path, bed_path)`

## Matrices: mean signal within each locus across replicates → mean within categories

For a replicate:
1. Open bigwig
2. Find intersections with L1HS bed file
3. Iterate through loci
    1. get mean signal at this location (populate dictionary for replicate)

Turn dictionary into df

Plot means across rows  (replicates)

## # Fold-change in full-length / truncated intergenic for each sample

For a bigwig:

1. Open file
2. Find intersections with L1HS intergenic bed file
3. Iterate through loci
    1. if full length add to full_length count
    2. if truncated add to truncate count
4. Fold-change is (full_length count/truncated count)
5. Store fold-change
