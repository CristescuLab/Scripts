#!/usr/bin/env bash

# This code transforms a fasta file to Cristescu-lab blast format and to dada2 format
#______________________________________
# Copyright (C) 2018  Jose Sergio Hleap
#______________________________________

# Arguments:
#   1. Fasta file
#   2. Path to accession2taxid lineages
#   3. Number of available cores
#   4. BOLD or NCBI
# The accessiion2taxid lineages is a file mapping the accession numbers, tax id,and lineages

set -e
#get the fasta prefix
acc2tax=$2
cpus=$3

process_NCBI()
{
# 1. fasta
# 2. acc2tax
# 3. cpu
#get accession list
prefix=${1%%.fasta}
seqkit -j $3 rmdup -s -i -m -o ${prefix}_processed.fasta -d duplicated.fa.gz \
-D duplicated.detail.txt $1
grep '>' ${prefix}_processed.fasta| cut -d' ' -f 1 | \
sed 's/>//' > ${prefix}.accessions
# get the lineages based on the accession numbers
LC_ALL=C awk 'FNR==NR{a[$1];next} ($1 in a)' ${prefix}.accesions $2 > ${prefix}.lineages
# get the mapping of accession number and stitle for renaming
format_them ${prefix}_processed.fasta ${prefix} $3
}


process_MIDORI()
{
# 1. zipped fasta
# 2. zipped taxon
# 3. cpus
prefix=${1%%.fasta.zip}
zcat $1| awk '{if ($0 ~ />/){ print $0 " MIDORI" ++count[$2]} else {print $0}}'| \
tee ${prefix}_numbered.fasta|  grep -a '>'| sed 's/ /\t/'|sed 's/>//'| \
awk -F $'\t' ' {OFS = FS} { t = $1; $1 = $2; $2 = t; print; } ' > ${prefix}.keyvalue
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' \
-k  ${prefix}.keyvalue ${prefix}_numbered.fasta| \
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' -k <(zcat $2 | sed 's/.__//g') | \
tee ${prefix}_lineages.fa | sed 's/[^ ]* />/' > ${prefix}_lineages4dada.fa
}


process_BOLD()
{
# 1. fasta
# 2. cpus
prefix=${1%%.fasta}
seqkit -j $1 grep -nrp COI $2| tee ${prefix}_COI_only.fasta| \
grep -a '>'| cut -d' ' -f 1 | sed 's/>//g'| > ${prefix}.accessions
python BOLD_lineage.py  ${prefix}.accessions ${prefix}.lineages $2
count=1
while [ $count -le 10 ]; do
mv failed.dump missing
nmiss=$(wc -l < missing)
chunk=$(( nmiss / 10 ))
if (( $chunk == 0 )); then
chunk=1
fi
python BOLD_lineage.py  missing ${prefix}.${count}.lineages $2 ${chunk}
(( count++ ))
done
#remove the not done
rm missing
sed 's/^/\^/' failed.dump | sed 's/$/ /' > missing.pattern
seqkit -j $2 grep -rnv -f missing.pattern $1 | \
seqkit -j $2 rmdup -s -i -m -o ${prefix}_processed.fasta -d ${prefix}_duplicated.fa.gz \
-D ${prefix}_duplicated.detail.txt
format_them ${prefix}_processed.fasta ${prefix} $2
}


format_them()
{
# 1. fasta
# 2. prefix
# 3. cpus
# 4. NCBI or BOLD
grep -a '>' $1 | sed 's/ /\t/'|sed 's/>//'| \
awk -F $'\t' ' {OFS = FS} { t = $1; $1 = $2; $2 = t; print; } ' > $2.keyvalue
if [[ "$4" == "NCBI" ]];then
  cut -f2,5 $2.lineages |sed 's/$/;/'| tr -s ';' > $2.map
else
  sed 's/$/;/' $2.lineages | tr -s ';' > $2.map
fi
# rename the fasta file
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' -k $2.keyvalue $1 | \
seqkit -j $3 replace -p ' (.+)$' -r ' {kv}' -k $2.map | tee $2_lineages.fa | \
sed 's/[^ ]* />/' > $2_lineages4dada.fa
}


if [[ "${@: -1}" == "NCBI" ]]; then
    process_NCBI $1 $2 $3
else
    process_BOLD $1 $2
fi