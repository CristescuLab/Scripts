#!/usr/bin/env bash
set -e
# 1. tsv blast file
# 2. query fasta
# 3. subject fasta
# 4. fasta out filename
count=0
while read line; do
#if [[ ${line} != *"qaccver"* ]]; then
  ((count++))
  read -a vars <<< $line
  echo "processing ${vars[0]} and ${vars[1]}"
  # get query sequence
  seqkit grep -p "${vars[0]}" $2 > tmp
  seqkit grep -p "${vars[1]}" $3 >> tmp
  mafft --quiet tmp | sed "s/>M0.*$/>Read${count}/"| sed 's/ /_/g' > tmp.aln
  trimal -nogaps -in tmp.aln >> tmp.fasta
#fi
done <$1
mafft --auto --thread -1 tmp.fasta | seqkit rmdup > $4
iqtree -pre ${4%%.fasta} -bb 1000 -alrt 0 -m TESTNEW  -con -s $4