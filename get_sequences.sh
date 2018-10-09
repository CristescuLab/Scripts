#!/usr/bin/env bash

# Usage bash get_sequences.sh <column of taxlevel> <name of filtered file>
## ALEXA
prefix=${2%%_filtered.tsv}
for i in `cut -f $1 $2 | sort -u| awk NF`;do
 ~/projects/def-mcristes/Software/seqkit grep -n -r -f <(grep "${i}" $2| \
 cut -f 1) ${prefix}.clean.fasta > ${prefix}_${i}.fasta
done

