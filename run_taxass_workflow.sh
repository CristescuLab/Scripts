#!/usr/bin/env bash
# To run this script:
# bash run_taxass_workflow.sh fastafile otutable prefix percentid database cpus
fasta=$1
otu_table=$2
prefix=$3
pid=$4
db=$5 #Greengenes or SILVAv132
cpus=$6
module load nixpkgs/16.09 gcc/7.3.0 boost/1.68.0 mothur/1.41.3 blast+/2.9.0
module load scipy-stack/2018b r-bundle-bioconductor/3.9
git clone https://github.com/McMahonLab/TaxAss
scripts=${PWD}/TaxAss/tax-scripts
cwd=$PWD

data2mothur(){
python3 - << EOF
import pandas as pd
df = pd.read_csv("${1}", sep='\t')
df = df.apply(lambda x: x / sum(x), axis=1).reset_index()
df.to_csv("${2}", sep='\t')
EOF
}

loop(){
# Step 5: filter BLAST results
Rscript ${scripts}/filter_seqIDs_by_pident.R ${prefix}_otus.custom.blast.table.modified \
${prefix}_ids.above.${1} ${1} TRUE
Rscript ${scripts}/filter_seqIDs_by_pident.R ${prefix}_otus.custom.blast.table.modified \
${prefix}_ids.below.${1} ${1} FALSE
# Step 7. recover missing seqIDs
python ${scripts}/find_seqIDs_blast_removed.py ${fasta} \
${prefix}_otus.custom.blast.table.modified ${prefix}_ids_${1}.missing
cat ${prefix}_ids.below.${1} ${prefix}_ids_${1}.missing > ${prefix}_ids.below.${1}.all
# Step 8. create fasta of seqIDs
python create_fastas_given_seqIDs.py ${prefix}_ids.above.${1} ${fasta} ${prefix}_otus.above.${1}.fasta
python create_fastas_given_seqIDs.py ${prefix}_ids.below.${1}.all ${fasta} ${prefix}_otus.below.${1}.fasta
# Step 9. assign taxonomy
mothur "#classify.seqs(fasta=${prefix}_otus.above.${1}.fasta, template=${ecosystem_fa},  taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=${cpus}, cutoff=0)"
mothur "#classify.seqs(fasta=${prefix}_otus.below.${1}.fasta, template=${general_fa}, taxonomy=${general_tax}, method=wang, probs=T, processors=${cpus}, cutoff=0)"
# Step 10. combine taxonomy files
cat ${prefix}_otus.above.${1}.custom.wang.taxonomy ${prefix}_otus.below.${1}.general.wang.taxonomy > ${prefix}_otus.${1}.taxonomy
# Step 12. reformat taxonomy files
sed 's/[[:blank:]]/\;/' <${prefix}_otus.${1}.taxonomy >${prefix}_otus.${1}.taxonomy.reformatted
mv ${prefix}_otus.${1}.taxonomy.reformatted ${prefix}_otus.${1}.taxonomy
sed 's/[[:blank:]]/\;/' <${prefix}_otus.general.taxonomy >${prefix}_otus.general.taxonomy.reformatted
mv ${prefix}_otus.general.taxonomy.reformatted ${prefix}_otus.general.taxonomy
# Step 13. compare taxonomy files
mkdir ${prefix}_conflicts_${1}
Rscript ${scripts}find_classification_disagreements.R ${prefix}_otus.${1}.taxonomy \
${prefix}_otus.general.taxonomy ${prefix}_ids.above.${1} conflicts_${1} ${1} 80 80
}


# adjust dadatest.R-generated otu table to mothur
sed '0,/^/{s/^/"name"\t/}' ${otu_table} > tmp
# Compute relative abundances
data2mothur tmp ${prefix}_otutab.tsv
# get FW files
unzip ./TaxAss/FreshTrain-files/FreshTrain18Aug2016.zip
# Set environmentalk variables to the fasta and taxonomy files
ecosystem_fa=FreshTrain18Aug2016/FreshTrain18Aug2016.fasta
ecosystem_tax=FreshTrain18Aug2016/FreshTrain18Aug2016.taxonomy
# Decompress the files for reference database
unzip ./TaxAss/FreshTrain-files/reshTrain-files/FreshTrain*${db}*
# Set the proper directory for the reference database
dirname=./TaxAss/FreshTrain-files/reshTrain-files/FreshTrain*${db}*
dirname=$(basename ${dirname})
# Get the prefix of the fileset
pre=${dirname%%.zip}
general_fa=${pre}/${pre}.fasta
general_tax=${pre}/${pre}.taxonomy
# Step 1.make BLAST database
makeblastdb -dbtype nucl -in ${ecosystem_fa} -input_type fasta -parse_seqids \
 -out custom.db
# Step 2. run BLAST
blastn -query ${fasta} -task megablast -db custom.db -out ${prefix}_otus.custom.blast \
-outfmt 11 -max_target_seqs 5
# Step 3. reformat BLAST results
blast_formatter -archive ${prefix}_otus.custom.blast \
-outfmt "6 qseqid pident length qlen qstart qend" -out ${prefix}_otus.custom.blast.table
# Step 4. recalculate pident
Rscript ${scripts}/calc_full_length_pident.R ${prefix}_otus.custom.blast.table \
${prefix}_otus.custom.blast.table.modified
# Step 5: filter BLAST results
Rscript ${scripts}/filter_seqIDs_by_pident.R ${prefix}_otus.custom.blast.table.modified \
${prefix}_ids.above.${pid} ${pid} TRUE
Rscript ${scripts}/filter_seqIDs_by_pident.R ${prefix}_otus.custom.blast.table.modified \
${prefix}_ids.below.${pid} ${pid} FALSE
# Step 6. check BLAST settings
Rscript ${scripts}/plot_blast_hit_stats.R ${prefix}_otus.custom.blast.table.modified \
${pid} plots
# Step 7. recover missing seqIDs
python ${scripts}/find_seqIDs_blast_removed.py ${fasta} \
${prefix}_otus.custom.blast.table.modified ${prefix}_ids.missing
cat ${prefix}_ids.below.${pid} ${prefix}_ids.missing > ${prefix}_ids.below.${pid}.all
# Step 8. create fasta of seqIDs
python create_fastas_given_seqIDs.py ${prefix}_ids.above.${pid} ${fasta} ${prefix}_otus.above.${pid}.fasta
python create_fastas_given_seqIDs.py ${prefix}_ids.below.${pid}.all ${fasta} ${prefix}_otus.below.${pid}.fasta
# Step 9. assign taxonomy
mothur "#classify.seqs(fasta=${prefix}_otus.above.${pid}.fasta, template=${ecosystem_fa},  taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=${cpus}, cutoff=0)"
mothur "#classify.seqs(fasta=${prefix}_otus.below.${pid}.fasta, template=${general_fa}, taxonomy=${general_tax}, method=wang, probs=T, processors=${cpus}, cutoff=0)"
# Step 10. combine taxonomy files
cat ${prefix}_otus.above.${pid}.custom.wang.taxonomy ${prefix}_otus.below.${pid}.general.wang.taxonomy > ${prefix}_otus.${pid}.taxonomy
# Step 11. general-only taxonomy
mothur "#classify.seqs(fasta=${fasta}, template=${general_fa}, taxonomy=${general_tax}, method=wang, probs=T, processors=${cpus}, cutoff=0)"
cp ${fasta%%.fasta}.general.wang.taxonomy > ${prefix}_otus.general.taxonomy
# Step 12. reformat taxonomy files
sed 's/[[:blank:]]/\;/' <${prefix}_otus.${pid}.taxonomy >${prefix}_otus.${pid}.taxonomy.reformatted
mv ${prefix}_otus.${pid}.taxonomy.reformatted ${prefiox}_otus.${pid}.taxonomy
sed 's/[[:blank:]]/\;/' <${prefix}_otus.general.taxonomy >${prefix}_otus.general.taxonomy.reformatted
mv ${prefix}_otus.general.taxonomy.reformatted ${prefix}_otus.general.taxonomy
# Step 13. compare taxonomy files
mkdir ${prefix}_conflicts_${pid}
Rscript ${scripts}/find_classification_disagreements.R ${prefix}_otus.${pid}.taxonomy \
${prefix}_otus.general.taxonomy ${prefix}_ids.above.${pid} ${prefix}_conflicts_${pid} ${pid} 80 80
# run multiple pident cutoffs
v=""
for i in seq $(( pid - 5 )) 100
do
  loop ${i}
  v=$(echo "${v} ${prefix}_conflicts_${i} ${prefix}_ids.above.${i} ${i}")
done
# Step 14. choose pident cutoff
Rscript ${scripts}/plot_classification_disagreements.R ${prefix}_otutab.tsv \
plots regular NA NA ${v}
# Step 15. make final taxonomy file
#Rscript ${scripts}/find_classification_disagreements.R ${prefix}_otus.${pid}.taxonomy
#otus.general.taxonomy ids.above.98 conflicts_${pid} ${pid} 80 80 final
#Rscript ${scripts}/find_classification_disagreements.R otus.${pid}.taxonomy \
#quickie ids.above.${pid} conflicts_${pid} ${pid} 80 80 final
#Rscript ${scripts}/plot_classification_disagreements.R ${prefix}_otutab.tsv \
# MakeSeqIDReadsOnly
#Rscript ${scripts}/find_classification_disagreements.R otus.${pid}.taxonomy \
#quickie ids.above.${pid} conflicts_${pid} ${pid} 80 80 final
## Step 15.5 plot TaxAss benefits
#Rscript ${scripts}/plot_classification_improvement.R final.taxonomy.pvalues \
#final.general.pvalues total.reads.per.seqID.csv plots final.taxonomy.names \
#final.general.names ids.above.${pid}
#
## Step 15.5.b. Plot improvement over using custom database alone
#mothur "#classify.seqs(fasta=${fasta}, template=${ecosystem_fa}, taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=32, cutoff=0)"
#cat otus.custom.wang.taxonomy > otus.custom.taxonomy
#sed 's/[[:blank:]]/\;/' <otus.custom.taxonomy >otus.custom.taxonomy.reformatted
#mv otus.custom.taxonomy.reformatted otus.custom.taxonomy
#mkdir conflicts_forcing
#Rscript ${scripts}/find_classification_disagreements.R otus.custom.taxonomy \
#otus.${pid}.80.80.taxonomy ids.above.${pid} conflicts_forcing NA 80 80 forcing
#Rscript ${scripts}/plot_classification_disagreements.R otus.abund plots \
#conflicts_forcing otus.custom.80.taxonomy otus.${pid}.80.80.taxonomy
#
## Step 16. tidy up
#bash ${scripts}/RunStep_16.sh

