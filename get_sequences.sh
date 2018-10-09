#!/usr/bin/env bash

# Usage bash get_sequences.sh <column of taxlevel> <name of filtered file>
## ALEXA
prefix=${2%%_filtered.tsv}
cut -f $1 $2 | sort -u| awk NF > ${prefix}.species
while read i;do
 ~/projects/def-mcristes/Software/seqkit grep -n -r -f <(grep "${i}" $2| \
 cut -f 1) ${prefix}.clean.fasta > ${prefix}_${i}.fasta
done<${prefix}.species

