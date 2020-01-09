#!/usr/bin/env bash
# script to subset a metagenome into metazoan bins
# Arguments:
# 1) db (including path)
# 2) fastq file (including path)
# 3) outprefix
# 4) path 2 pygrep
# 5) cpus

# exit when any command fails
set -e
find_euks(){
python3 - << EOF
import pandas as pd
import ete3

ncbi = ete3.NCBITaxa()
levels = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus',
          'species']
def get_lineage(taxid):
    taxid_lineage = {v: k for k, v in ncbi.get_rank(
        ncbi.get_lineage(taxid)).items() if v != 'no rank'}
    lineage = ncbi.translate_to_names(
        [taxid_lineage[x] for x in levels if x in taxid_lineage])

    return 'Eukaryota' in lineage

cols = ['class', 'qid', 'krakid', 'length', 'kmerid']
df = pd.read_csv("$1", sep='\t', header=None, names=cols)
df = df[df.krakid != 0]
df['Euk'] = df.krakid.apply(get_lineage)
euks = df[df.Euk]
euks.qid.to_csv("$2", header=False, index=False)
EOF
}

dbname=$(basename $1)
prefix="${3}"_"${dbname}"
# Run Kaken classification
kraken2 --db "${1}" <(seqkit fq2fa "${2}") > "${prefix}".kraken2
# Get the unclassified
grep "^U" "${prefix}".kraken2| cut -f 2 > "${prefix}".unclassified
python3 ${4} ${2} -c ${5} -p "${prefix}".unclassified > "${prefix}".unclassified.fq
# get the eukaryotes from the classified sequences
find_euks "${prefix}".kraken2 "${prefix}".euks
python3 ${4} ${2} -c ${5} -p "${prefix}".euks > "${prefix}".euks.fq

