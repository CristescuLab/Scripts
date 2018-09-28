#!/usr/bin/env bash
#TODO: download only if necessary. check MDsums and avoid duplications in the final file
tmp_folder=$1
cpus=4
echo bye | lftp -e "mirror -c" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
echo "will cite" | parallel --citation
for i in *.gz; do
    a=`zcat ${i} | cut -f1-3 | tee tmp| tail -n 1| cut -f 1`
    if LC_ALL=C fgrep -q -m 1 blah accession2taxid; then
    echo "$i has previously been done"; else
    parallel -j ${cpus} --tmpdir ${tmp_folder} --pipepart -a tmp --block 500m \
    'taxonkit lineage -i 3 | taxonkit reformat -i 4' >> accession2taxid
    rm tmp
    fi
done
sort -T ${tmp_folder} -u accession2taxid > accession2taxid.lineage
