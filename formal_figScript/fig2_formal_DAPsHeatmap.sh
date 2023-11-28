#!/bin/bash

#SBATCH --job-name=deeptools
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=6hours

ml purge
ml deepTools

# for panel H

#bamCoverage -b /scicore/home/tschoppp/wang0007/scATAC/peakcall/somite/somite_clssub_4.bam -o /scicore/home/tschoppp/wang0007/scATAC/peakcall/somite/somite_clssub_4.cpm.bw --effectiveGenomeSize 1060000000 --normalizeUsing CPM

#bamCoverage -b /scicore/home/tschoppp/wang0007/scATAC/peakcall/nasal/nasal_clssub_6.bam -o /scicore/home/tschoppp/wang0007/scATAC/peakcall/nasal/nasal_clssub_6.cpm.bw --effectiveGenomeSize 1060000000 --normalizeUsing CPM

#bamCoverage -b /scicore/home/tschoppp/wang0007/scATAC/peakcall/limb/limb_clssub_5.bam -o /scicore/home/tschoppp/wang0007/scATAC/peakcall/limb/limb_clssub_5.cpm.bw --effectiveGenomeSize 1060000000 --normalizeUsing CPM

#computeMatrix scale-regions -S /scicore/home/tschoppp/wang0007/scATAC/peakcall/somite/somite_clssub_4.cpm.bw \
#                                 /scicore/home/tschoppp/wang0007/scATAC/peakcall/nasal/nasal_clssub_6.cpm.bw  \
#                                 /scicore/home/tschoppp/wang0007/scATAC/peakcall/limb/limb_clssub_5.cpm.bw \
#                              -R daps_chondrocytes_Top500limb.bed daps_chondrocytes_Top500nasal.bed daps_chondrocytes_Top500somite.bed \
#                              --beforeRegionStartLength 500 \
#                              --afterRegionStartLength 500 \
#                              -o daps_chondrocytes_Top500.cpm.mat.gz

#plotHeatmap --colorMap jet -m daps_chondrocytes_Top500.cpm.mat.gz \
#      -out daps_chondrocytes_Top500_changeColor.cpm.pdf

plotHeatmap --colorMap YlOrBr -m daps_chondrocytes_Top500.cpm.mat.gz \
      --xAxisLabel "peak distance (bp)" --startLabel "start" --endLabel "end" \
      -out daps_chondrocytes_Top500_colYlOrBr.cpm.pdf


