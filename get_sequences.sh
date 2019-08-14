#!/usr/bin/env bash
# Need to load muscle:
# module load nixpkgs/16.09  gcc/5.4.0 mafft/7.310
# Usage bash get_sequences.sh <column of taxlevel> <filtered file> <fasta file> <cpus>
column=$1
filtered_file=$2
fas=$3
cpus=$4

prefix=${filtered_file%%_filtered.tsv}
cut -f ${column} ${filtered_file} | sort -u| awk NF > ${prefix}.species
while read i;do
sp=`echo ${i}| tr ' ' '_'`
 outfn=`echo ${prefix}_${i}.fas|tr ' ' '_'`
 if [[ ! -f ${fas%%.fasta}_${sp}.fasta ]]; then
    seqkit -j ${cpus} grep -n -r -f <(grep "${i}" ${filtered_file}| cut -f 1) \
    ${fas} > ${prefix}_${sp}.fasta
 fi
mafft --auto --quiet --thread ${cpus} ${prefix}_${sp}.fasta > ${prefix}_${sp}
done<${prefix}.species

