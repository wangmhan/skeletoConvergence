#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --qos=1day

ml purge #ml avail
ml load SAMtools #1.3
ml load MACS2 #2.2.5

while getopts ":b:o:" options; do # A : after the variable means that we are expecting a value for it
        case $options in
        b ) bamfile=$OPTARG;;
		o ) output=$OPTARG;;
        \? ) echo "Invalid option -$OPTARG"
        exit 1;;
        :) echo "Option -$OPTARG requires an argument"
        exit 1;;
        esac
done

samtools index $bamfile
macs2 callpeak -t $bamfile -f BAMPE -n $output -g 1.06e9 --keep-dup all --nomodel --shift 100 --extsize 200 --call-summits



