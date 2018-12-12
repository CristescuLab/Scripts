#!/usr/bin/env bash


cpus=32
# check if the primers still there in any configuration
seqkit -j 32 locate -d -p "${mlCOIintFprimer}" ${outdir}/${data}_pear.assembled.fastq > fwd
seqkit -j 32 locate -d -p "${mlCOIintFprimer}" ${outdir}/${data}_pear.assembled.fastq > rev

if [[ $(wc -l <fwd) -ge 2 || $(wc -l <rev) -ge 2 ]]; then
    cutadapt -j ${cpus} -a ${Adapter2rc} -g ${mlCOIintFprimer} \
    -g ${Adapter2rc} -a ${mlCOIintFprimer} -m 100 -n 2 --match-read-wildcards \
    ${outdir}/${data}_pear.assembled.fastq > ${outdir}/${data}_pear.assembled_secondcut.fastq
fi
