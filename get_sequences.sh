#!/usr/bin/env bash

# Usage bash get_sequences.sh <column of taxlevel> <name of filtered file>
prefix=${2%%_filtered.tsv}
cut -f $1 $2 | sort -u| awk NF > ${prefix}.species
while read i;do
 outfn=`echo ${prefix}_${i}.fas|tr ' ' '_'`
 ~/projects/def-mcristes/Software/seqkit grep -n -r -f <(grep "${i}" $2| \
 cut -f 1) ${prefix}.clean.fasta > ${outfn}
done<${prefix}.species

