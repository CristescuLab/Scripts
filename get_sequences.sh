#!/usr/bin/env bash
# Need to load muscle:
# module load nixpkgs/16.09  gcc/5.4.0 intel/2016.4 muscle/3.8.31
# Usage bash get_sequences.sh <column of taxlevel> <name of filtered file>
prefix=${2%%_filtered.tsv}
cut -f $1 $2 | sort -u| awk NF > ${prefix}.species
while read i;do
 outfn=`echo ${prefix}_${i}.fas|tr ' ' '_'`
 seqkit grep -n -r -f <(grep "${i}" $2| cut -f 1) ${prefix}.clean.fasta | \
 muscle> ${outfn}
done<${prefix}.species

