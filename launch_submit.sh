#!/usr/bin/env bash
# This script is a simple for loop to launch multiple scripts
# jobs. It requires that the fasta to be blast be in a subfolder within where you will
# be launching the job. It assumes that your folder is named with the last dot separated
# field being of your sample. The submit_blast_n_process.sh
# should be in your folder or in your path
# Usage:
# bash launch_submit.sh <pattern_folder> <database_with_path> <path_to_code> <cpus> <mem> <account>
# <cpus> <mem> <account> are optional, defaulting to 8, 32G and def-mcristes
pattern_folder=$1
db=$2
path_to_code=$3
cpus=$4
mem=$5
acc=$6


if [ -z ${cpus} ]; then cpus=8; fi
if [ -z ${mem} ]; then mem=32G; fi
if [ -z ${acc} ]; then acc=def-mcristes; fi

for i in `find . -type d -name "*${pattern_folder}*"`; do
name=`echo ${i} | rev | cut -d'.' -f 1 | rev`
sbatch --job-name=${name} -o ${name}.out --cpus-per-task=${cpus}\
--export=file_path=${i},db_path=${db},cpus=${cpus} --mem=${mem} \
-A ${acc} ${path_to_code}
done