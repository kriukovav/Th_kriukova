#!/bin/bash

#SBATCH --cpus-per-task=40         # Run on 40 CPU
#SBATCH --mem=16gb                 # Job memory request
#SBATCH --partition=medium

# read a line with sample names
# $line - a line containing path to R1 in the first column, path to R2 in the second column, and sample_id in the third column

R1=$(echo "$i" | cut -f 1)
R2=$(echo "$i" | cut -f 2)
sample_id=$(echo "$i" | cut -f 3)
    
#body of script. mixcr analyse amplicon.
echo $R1
echo $R2
echo $sample_id
$mixcr_path analyze milab-human-tcr-rna-race-cdr3 \
-Massemble.consensusAssemblerParameters.assembler.minRecordsPerConsensus=2 \
-MrefineTagsAndSort.parameters.postFilter=null \
$R1 $R2 $sample_id 

echo "${sample_id} finished" > ./log_finished/${sample_id}_finished.txt
