#!/usr/bin/env bash
# Need to load muscle:
# module load nixpkgs/16.09  gcc/5.4.0 mafft/7.310
# Usage bash get_sequences.sh <column of taxlevel> <name of filtered file> <cpus>
prefix=${2%%_filtered.tsv}
cut -f $1 $2 | sort -u| awk NF > ${prefix}.species
while read i;do
 outfn=`echo ${prefix}_${i}.fas|tr ' ' '_'`
 if [ ! -f ${prefix}_sub.fasta ]; then
    seqkit -j $3 grep -n -r -f <(grep "${i}" $2| cut -f 1) \
    ${prefix}.clean.fasta > ${prefix}_sub.fasta
 fi
mafft --auto --quiet --thread $3 ${prefix}_sub.fasta > ${outfn}
rm ${prefix}_sub.fasta
done<${prefix}.species

