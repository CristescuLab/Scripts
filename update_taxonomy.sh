#!/usr/bin/env bash
#TODO: download only if necessary. check MDsums and avoid duplications in the final file
cpus=4
echo bye | lftp -e "mirror -c" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
if [[ -f tmp ]]; then
rm tmp
fi
echo "will cite" | parallel --citation
for i in *.gz; do
    zcat ${i} | cut -f1-3 > tmp
    parallel -j ${cpus} --pipepart -a tmp --block 500m \
    'taxonkit lineage -i 3 | taxonkit reformat -i 4' >> accession2taxid
    rm tmp
done
sort -u accession2taxid > accession2taxid.lineage
#taxonkit -j `nproc` lineage -i 3 accession2taxid | \
#taxonkit -j `nproc` reformat -i 4 > accession2taxid.lineage