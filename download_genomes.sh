#!/usr/bin/env bash
domain=$1

download_one(){
ftp_path=$1
prefix=$(basename ${ftp_path})
wget -c ${ftp_path}/${prefix}_genomic.fna.gz
wget -c ${ftp_path}/md5checksums.txt -O ${prefix}_md5checksums.txt
md5sum -c <(grep GCF_000002515.2_ASM251v1_genomic.fna.gz ${prefix}_md5checksums.txt)
}

tid=$(echo ${domain}| taxonkit name2taxid| cut -f 2)
taxonkit list --ids 2759 --indent "" > ${domain}_species_taxid
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
sed 1d assembly_summary_refseq.txt | sed -e '1s/^# //' -e 's/"/$/g' > assembly_summary.tsv
csvtk -t pretty assembly_summary.tsv > assembly_summary.tsv.pretty
csvtk cut -t -f organism_name assembly_summary.tsv| cut -d ' ' -f 1,2 | \
csvtk freq -t -n -r | csvtk pretty -t > assembly.sps
csvtk grep -t -f species_taxid -P ${domain}_species_taxid assembly_summary.tsv \
| csvtk grep -t -f assembly_level -i -p "Complete Genome" > ${domain}_bytaxid.tsv

csvtk cut -t -f ftp_path ${domain}_bytaxid.tsv| sed 1d | \
rush -v prefix='{}/{%}' \
        ' \
            wget -c {prefix}_genomic.fna.gz; \
            wget -c {prefix}_genomic.gbff.gz; \
            wget -c {prefix}_genomic.gff.gz; \
            wget -c {prefix}_cds_from_genomic.fna.gz \
        ' \
        -j 10 -c -C download.rush