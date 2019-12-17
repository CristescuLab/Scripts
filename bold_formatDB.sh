#!/usr/bin/env bash
# 1. number of cpus
# 2. Fasta file from BOLD

prefix=${2%%.fasta}

seqkit -j $1 grep -nrp COI $2| tee ${prefix}_COI_only.fasta| \
grep -a '>'| sed 's/ /\t/'| cut -f 2| sed 's/ sp.*$//'  | \
cut -d'|' -f 1| sed 's/_incertae_sedis//'| sed 's/\S*\(-SUPPRESSED\)\S*//g'| \
sed 's/ .*[0-9].*//'| sed 's/GEN.*$//'| sed 's/group$//'| sed 's/_grp//' | \
sed 's/.(nomen nudum)//'| sed 's/^.*(//'| sed 's/)//'|sed 's/NZD//'| \
sed 's/_gen.*$//'| sed 's/TW//'| sed 's/_order//' | \
taxonkit -j $1 name2taxid| tee names2taxid | awk ' $2 == "" '| \
sort -u >${prefix}.missing 

taxonkit -j $1 lineage -i 2 names2taxid | taxonkit -j $1 reformat -i 3

LC_ALL=C awk 'FNR==NR{a[$1];next} ($2 in a)' ${prefix}.missing ${prefix}_COI_only.fasta| \
cut -d' ' -f1 | sed 's/>//'| sort -u > missingBOLDacc.list

python BOLD_lineage.py missingBOLDacc.list missingBOLDlineages.list 10
comm -13 <(cut -f1 missingBOLDlineages.list|sort) missingBOLDacc.list > still_missing
count=1
while [[ `wc -l < still_missing` -ne 0 ]]; do
python BOLD_lineage.py still_missing tmp 10
if [ `wc -l < tmp` -ne 0 ]; then
cat tmp >> missingBOLDlineages.list
fi
comm -13 <(cut -f1 missingBOLDlineages.list|sort) missingBOLDacc.list > still_missing
count=$(( count + 1 ))
done



grep -a '>' BOLD_Metazoa_split_COI_only.fasta| cut -d' ' -f1| sed 's/>//'| sort -u > BOLD_ACC.list
