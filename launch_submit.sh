#!/usr/bin/env bash
# This script is a simple for loop to launch multiple submit_blast_n_process.sh
# jobs. It requires that the fasta to be blast be in a subfolder within where you will
# be launching the job. It assumes that your folder is names into at least 5 field,
# separated by ., where the 5th field if the name of your sample. The submit_blast_n_process.sh
# should be in your folder or in your path
# Usage:
# bash launch_submit.sh <pattern_folder> <database_with_path>

for i in `find . -type d -name *"$1"*`; do
name=`cut -d'.' -f5 ${i}`
sbatch --job-name=${name} --export=file_path=${i},db_path=$2 submit_blast_n_process.sh
done