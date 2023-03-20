#!/bin/bash

#SBATCH --job-name=subset-bam
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=6hours

ml purge

while getopts ":b:c:o:" options; do # A : after the variable means that we are expecting a value for it
        case $options in
        b ) bamfile=$OPTARG;;
		c ) cellid=$OPTARG;;
		o ) output=$OPTARG;;
        \? ) echo "Invalid option -$OPTARG"
        exit 1;;
        :) echo "Option -$OPTARG requires an argument"
        exit 1;;
        esac
done

/scicore/home/tschoppp/wang0007/software/subset-bam/subset-bam \
	--bam $bamfile \
	--cell-barcodes $cellid \
	--out-bam $output

# sbatch macs2_peakcall.sh -b ~/scATAC/reslt_scATAC/forBlockcourse/rawdata/possorted_bam.bam -c ~/scATAC/reslt_scATAC/forBlockcourse/results/cellid_neuralCrest.txt -o cellid_neuralCrest.bam



