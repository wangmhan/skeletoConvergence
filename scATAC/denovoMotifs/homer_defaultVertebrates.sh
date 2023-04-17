#!/bin/bash

#SBATCH --job-name=homer
#SBATCH --time=3-0:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --qos=1week

ml purge

while getopts ":b:o:" options; do # A : after the variable means that we are expecting a value for it
        case $options in
        	b ) bedfile=$OPTARG;;
		o ) output=$OPTARG;;
        \? ) echo "Invalid option -$OPTARG"
        exit 1;;
        :) echo "Option -$OPTARG requires an argument"
        exit 1;;
        esac
done

findMotifsGenome.pl $bedfile \
	/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
	$output -len 8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
	-mset vertebrates -size -250,250  -fdr 5 -p 4 -cache 4000



