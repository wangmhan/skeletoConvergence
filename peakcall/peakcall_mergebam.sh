#!/bin/bash

#SBATCH --job-name=merge-bam
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=6hours

ml purge
ml load SAMtools

while getopts ":o:m:n:" options; do # A : after the variable means that we are expecting a value for it
        case $options in
        o ) output=$OPTARG;;
	m ) bamfile1=$OPTARG;;
	n ) bamfile2=$OPTARG;;
	\? ) echo "Invalid option -$OPTARG"
        exit 1;;
        :) echo "Option -$OPTARG requires an argument"
        exit 1;;
        esac
done

samtools merge $output $bamfile1 $bamfile2
#echo "$output"
#echo "$bamfile1"

# sbatch merge_bam.sh -o merged.bam -m bam1.bam -n bam2.bam




