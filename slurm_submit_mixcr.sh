#!/bin/bash

# this script executes mixcr analyze script in slurm batches

mixcr_path=$1 #please specify path to mixcr software when you execute this script

mkdir -p ./outs/mixcr/log_finished
mkdir -p ./outs/mixcr/stdout_files
cd ./outs/mixcr

while read i ; do
    sample_id=$(echo "$i" | cut -f 3)
    sbatch -J $sample_id --export=i="$i",mixcr_path="$mixcr_path" --output ./stdout_files/$sample_id.txt ../../run_mixcr.sh 
done < ../../filenames.txt

number_submited=$(cat ../../filenames.txt | grep -c ^)
number_finished=$(ls ./log_finished | grep -c .txt)

while [ $number_submited -gt $number_finished ]; do
  sleep 60
  number_submited=$(cat ../../filenames.txt | grep -c ^)
  number_finished=$(ls ./log_finished | grep -c  finished.txt)
done

rm ./log_finished/*