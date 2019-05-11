#!/usr/bin/env bash

tmp_folder=$1
cpus=4
tmp=${tmp_folder}/tmp
echo bye | lftp -e "mirror -c" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
echo "will cite" | parallel --citation
for i in *.gz; do
    #a=`zcat ${i} | tail -n +2| cut -f1-3 | tee ${tmp}| tail -n 1| cut -f1,2`
    a=`zcat ${i} | tail -1| cut -f1,2`
    if LC_ALL=C grep -q -m 1 "^${a}" accession2taxid.lineage; then
    echo "$i has previously been done"; else
    echo "Processing ${i}"
    zcat ${i} | tail -n +2| cut -f1-3 > ${tmp}
    parallel -j ${cpus} --tmpdir ${tmp_folder} --pipepart \
    -a ${tmp} --block 500m \
    'taxonkit lineage -i 3 | taxonkit reformat -i 4' >> accession2taxid.lineage
    rm ${tmp}
    fi
done
sort -T ${tmp_folder} -u accession2taxid.lineage > ${tmp} && \
mv ${tmp} accession2taxid.lineage && rm ${tmp} tmp
