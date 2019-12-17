#!/usr/bin/env bash

fasta=$1
otu_table=$2
prefix=$3
pid=$4
db=$5 #Greengenes or SILVAv132
module load nixpkgs/16.09 gcc/7.3.0 boost/1.68.0 mothur/1.41.3 blast+/2.9.0
module load scipy-stack/2018b
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
Rscript ${scripts}/filter_seqIDs_by_pident.R otus.custom.blast.table.modified \
ids.above.${1} ${1} TRUE
Rscript ${scripts}/filter_seqIDs_by_pident.R otus.custom.blast.table.modified \
ids.below.${1} ${1} FALSE
# Step 7. recover missing seqIDs
python ${scripts}/find_seqIDs_blast_removed.py otus.fasta \
otus.custom.blast.table.modified ids.missing
cat ids.below.${1} ids.missing > ids.below.${1}.all
# Step 8. create fasta of seqIDs
python create_fastas_given_seqIDs.py ids.above.98 otus.fasta otus.above.${1}.fasta
python create_fastas_given_seqIDs.py ids.below.98.all otus.fasta otus.below.${1}.fasta
# Step 9. assign taxonomy
mothur "#classify.seqs(fasta=otus.above.${1}.fasta, template=${ecosystem_fa},  taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=32, cutoff=0)"
mothur "#classify.seqs(fasta=otus.below.${1}.fasta, template=${general_fa}, taxonomy=${general_tax}, method=wang, probs=T, processors=32, cutoff=0)"
# Step 10. combine taxonomy files
cat otus.above.${1}.custom.wang.taxonomy otus.below.${1}.general.wang.taxonomy > otus.${1}.taxonomy
# Step 12. reformat taxonomy files
sed 's/[[:blank:]]/\;/' <otus.${1}.taxonomy >otus.${1}.taxonomy.reformatted
mv otus.${1}.taxonomy.reformatted otus.${1}.taxonomy
sed 's/[[:blank:]]/\;/' <otus.general.taxonomy >otus.general.taxonomy.reformatted
mv otus.general.taxonomy.reformatted otus.general.taxonomy
# Step 13. compare taxonomy files
mkdir conflicts_${1}
Rscript ${scripts}find_classification_disagreements.R otus.${1}.taxonomy \
otus.general.taxonomy ids.above.${1} conflicts_${1} ${1} 80 80
}


# adjust dadatest.R-generated otu table to mothur
sed '0,/^/{s/^/"name"\t/}' ${otu_table} > tmp
data2mothur tmp ${prefix}_otutab.tsv
# get FW files
unzip ./TaxAss/FreshTrain-files/FreshTrain18Aug2016.zip
ecosystem_fa=FreshTrain18Aug2016/FreshTrain18Aug2016.fasta
ecosystem_tax=FreshTrain18Aug2016.taxonomy
unzip ./TaxAss/FreshTrain-files/reshTrain-files/FreshTrain*${db}*
dire=./TaxAss/FreshTrain-files/reshTrain-files/FreshTrain*${db}*
dire=$(basename ${dirname})
pre=${dire%%.zip}
general_fa=${pre}/${pre}.fasta
general_tax=${pre}/${pre}.taxonomy
# Step 1.make BLAST database
makeblastdb -dbtype nucl -in ${ecosystem_fa} -input_type fasta -parse_seqids \
 -out custom.db
# Step 2. run BLAST
blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast \
-outfmt 11 -max_target_seqs 5
# Step 3. reformat BLAST results
blast_formatter -archive otus.custom.blast \
-outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table
# Step 4. recalculate pident
Rscript ${scripts}/calc_full_length_pident.R otus.custom.blast.table \
otus.custom.blast.table.modified
# Step 5: filter BLAST results
Rscript ${scripts}/filter_seqIDs_by_pident.R otus.custom.blast.table.modified \
ids.above.${pid} ${pid} TRUE
Rscript ${scripts}/filter_seqIDs_by_pident.R otus.custom.blast.table.modified \
ids.below.${pid} ${pid} FALSE
# Step 6. check BLAST settings
Rscript ${scripts}/plot_blast_hit_stats.R otus.custom.blast.table.modified \
${pid} plots
# Step 7. recover missing seqIDs
python ${scripts}/find_seqIDs_blast_removed.py otus.fasta \
otus.custom.blast.table.modified ids.missing
cat ids.below.${pid} ids.missing > ids.below.${pid}.all
# Step 8. create fasta of seqIDs
python create_fastas_given_seqIDs.py ids.above.98 otus.fasta otus.above.${pid}.fasta
python create_fastas_given_seqIDs.py ids.below.98.all otus.fasta otus.below.${pid}.fasta
# Step 9. assign taxonomy
mothur "#classify.seqs(fasta=otus.above.${pid}.fasta, template=${ecosystem_fa},  taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=32, cutoff=0)"
mothur "#classify.seqs(fasta=otus.below.${pid}.fasta, template=${general_fa}, taxonomy=${general_tax}, method=wang, probs=T, processors=32, cutoff=0)"
# Step 10. combine taxonomy files
cat otus.above.${pid}.custom.wang.taxonomy otus.below.${pid}.general.wang.taxonomy > otus.${pid}.taxonomy
# Step 11. general-only taxonomy
cp otus.below.${pid}.general.wang.taxonomy otus.general.taxonomy
# Step 12. reformat taxonomy files
sed 's/[[:blank:]]/\;/' <otus.${pid}.taxonomy >otus.${pid}.taxonomy.reformatted
mv otus.${pid}.taxonomy.reformatted otus.${pid}.taxonomy
sed 's/[[:blank:]]/\;/' <otus.general.taxonomy >otus.general.taxonomy.reformatted
mv otus.general.taxonomy.reformatted otus.general.taxonomy
# Step 13. compare taxonomy files
mkdir conflicts_${pid}
Rscript ${scripts}/find_classification_disagreements.R otus.${pid}.taxonomy \
otus.general.taxonomy ids.above.${pid} conflicts_${pid} ${pid} 80 80
# run multiple pident cutoffs
for i in 90 92 94 96 98 99 100
do
  loop ${i}
done
# Step 14. choose pident cutoff
Rscript ${scripts}/plot_classification_disagreements.R ${prefix}_otutab.tsv \
plots regular NA NA conflicts_90 ids.above.90 90 conflicts_92 ids.above.92 92 \
conflicts_94 ids.above.94 94 conflicts_96 ids.above.96 96 conflicts_98 \
ids.above.98 98 conflicts_99 ids.above.99 100 conflicts_100 ids.above.100 100
# Step 15. make final taxonomy file
Rscript ${scripts}/find_classification_disagreements.R otus.98.taxonomy
otus.general.taxonomy ids.above.98 conflicts_${pid} ${pid} 80 80 final
Rscript ${scripts}/find_classification_disagreements.R otus.${pid}.taxonomy \
quickie ids.above.${pid} conflicts_${pid} ${pid} 80 80 final
Rscript ${scripts}/plot_classification_disagreements.R ${prefix}_otutab.tsv \
 MakeSeqIDReadsOnly
Rscript ${scripts}/find_classification_disagreements.R otus.${pid}.taxonomy \
quickie ids.above.${pid} conflicts_${pid} ${pid} 80 80 final
# Step 15.5 plot TaxAss benefits
Rscript ${scripts}/plot_classification_improvement.R final.taxonomy.pvalues \
final.general.pvalues total.reads.per.seqID.csv plots final.taxonomy.names \
final.general.names ids.above.${pid}

# Step 15.5.b. Plot improvement over using custom database alone
mothur "#classify.seqs(fasta=otus.fasta, template=${ecosystem_fa}, taxonomy=${ecosystem_tax}, method=wang, probs=T, processors=32, cutoff=0)"
cat otus.custom.wang.taxonomy > otus.custom.taxonomy
sed 's/[[:blank:]]/\;/' <otus.custom.taxonomy >otus.custom.taxonomy.reformatted
mv otus.custom.taxonomy.reformatted otus.custom.taxonomy
mkdir conflicts_forcing
Rscript ${scripts}/find_classification_disagreements.R otus.custom.taxonomy \
otus.${pid}.80.80.taxonomy ids.above.${pid} conflicts_forcing NA 80 80 forcing
Rscript ${scripts}/plot_classification_disagreements.R otus.abund plots \
conflicts_forcing otus.custom.80.taxonomy otus.${pid}.80.80.taxonomy

# Step 16. tidy up
bash ${scripts}/RunStep_16.sh

