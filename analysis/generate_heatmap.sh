#!/bin/bash



input_file=$1
#prefix is base name of input file of input
prefix=$(basename $input_file)

right_now=$(date +"%H:%M:%S")
echo "***Processing $input_file files at $right_now***"


trunc_inter="data/truncate_intergenic.bed"
trunc_intra="data/truncate_intronic.bed"
full_inter="data/full_intergenic.bed"
full_intra="data/full_intronic.bed"

computeMatrix reference-point \
--regionsFileName $full_inter $trunc_inter $full_intra $trunc_intra   \
--scoreFileName $input_file \
--referencePoint TES \
--binSize 10 \
--upstream 7000 \
--downstream 5000 \
--sortUsing region_length \
--outFileName heatmaps/${prefix}.gz


echo "computeMatrix completed: files at heatmaps/${prefix}.gz"

plotHeatmap -m heatmaps/${prefix}.gz \
--outFileName heatmaps/${prefix}.png \
--missingDataColor white \
--colorList 'purple,blue' \
--colorNumber 256 \
--sortUsing region_length \
--whatToShow 'plot, heatmap and colorbar' \
--heatmapWidth 7 --heatmapHeight 30

echo "plotHeatmap completed: files at heatmaps/${prefix}.png"
