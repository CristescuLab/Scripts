#!/usr/bin/env bash
# Arguments
# 1) fileset (prefix of R1 and R2)
# 2) preferred prefix
# 3) database
# 4) usearch executable
# 5) Forward primer
# 6) Reverse primer
# 7) primer name
# 8) Number of threads for blast (ONLY ONE IF USING GNU PARALLEL!!)

length_stats(){
python - << EOF
from collections import Counter
lengths = {}
name = None
seq = ''
with open("$1") as F:
    for line in F:
        if line.startswith('>'):
            if name is not None:
                lengths[name] = len(seq.strip().replace('\n',''))
            seq = ''
            name = line.strip()
        else:
            seq += line
    if name not in lengths:
        lengths[name] = len(seq.strip().replace('\n',''))
c = Counter(lengths.values())
mc = c.most_common(1)[0][0]
gt = sum([x[1] for x in c.items() if x[0] > mc ])
lt = sum([x[1] for x in c.items() if x[0] < mc ])
line = '\n%d sequences are longer than the mode\n%d sequences are ' \
       'shorter than the mode' % (gt, lt)
sor = sorted(c.items(), key=lambda x: x[0])
with open('${2}_distribution.txt', 'w') as F:
    F.write('\n'.join(['%d:%d' % x for x in sor]))
    F.write(line)
EOF
}

blast(){
# arguments
# 1) File
# Out prefix (with path)
blastn -db ${db} -query $1 -evalue 0.0001 -perc_identity 90 -outfmt \
"6 qseqid sseqid pident evalue qcovs qlen length staxid stitle" \
-max_target_seqs 40 -num_threads ${3} > ${2}_nt.hits && echo -e "\n#END" >> ${2}_nt.hits
}
# Define variables
fileset=$1
prefix=$2
#`echo ${fileset}|rev|cut -d'.' -f1|rev`
outdir=${prefix}_output
db=$3
u=$4
#COI: GGWACWGGWTGAACWGTWTAYCCYC TAAACTTCAGGGTGACCAAAAAATCA
#12S: GTCGGTAAAACTCGTGCCAGC CATAGTGGGGTATCTAATCCCAGTTTG
ADAPTER_FWD=$5 #"GGWACWGGWTGAACWGTWTAYCCYCC"
ADAPTER_REV=$6 #"TAAACTTCAGGGTGACCAAAAAATCA"
primer_name=_$7
Adapter1rc=`echo $ADAPTER_FWD | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
Adapter2rc=`echo $ADAPTER_REV | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
# create output folder
mkdir -p ${outdir}${primer_name}
outdir=${outdir}${primer_name}
# run first cutadapt
cutadapt -g "${ADAPTER_FWD}" -G "${ADAPTER_REV}" -a "${Adapter2rc}" \
-A "${Adapter1rc}" -o ${outdir}/${prefix}.1.fastq.gz \
-p ${outdir}/${prefix}.2.fastq.gz -m 130 --match-read-wildcards -q 25 --trim-n \
-n 2 --untrimmed-output ${outdir}/untrimmed.${prefix}.fastq \
--untrimmed-paired-output ${outdir}/untrimmed.paired.${prefix}.fastq \
${fileset}_R1.fastq.gz ${fileset}_R2.fastq.gz > ${outdir}/${prefix}.log1
# Merge reads
pear -f ${outdir}/${prefix}.1.fastq.gz -r  ${outdir}/${prefix}.2.fastq.gz \
-o ${outdir}/${prefix}_pear -q 25 -t 100 -s 2 > ${outdir}/${prefix}_pear.log
# check if adapters still there
seqkit -j 8 locate -d -p "${ADAPTER_FWD}" \
${outdir}/${prefix}_pear.assembled.fastq > ${outdir}/${prefix}.fwd
seqkit -j 8 locate -d -p "${ADAPTER_REV}" \
${outdir}/${prefix}_pear.assembled.fastq > ${outdir}/${prefix}.rev
# cut the adapter if they remain
if [[ $(wc -l <${outdir}/${prefix}.fwd) -ge 2 || $(wc -l <${outdir}/${prefix}.rev) -ge 2 ]]; then
cutadapt -a $Adapter2rc -m 100 -n 2 --too-short-output ${outdir}/${prefix}.short.fastq \
-o ${outdir}/${prefix}.3trimmed.fastq ${outdir}/${prefix}_pear.assembled.fastq > ${outdir}/${prefix}.log2
cutadapt -g ${ADAPTER_FWD} -n 2 -o ${outdir}/${prefix}.5trimmed.fastq \
${outdir}/${prefix}.3trimmed.fastq > ${outdir}/${prefix}.log3; else
cp ${outdir}/${prefix}_pear.assembled.fastq ${outdir}/${prefix}.3trimmed.fastq
fi
# dereplicate
$u -fastx_uniques ${outdir}/${prefix}.3trimmed.fastq -fastaout \
${outdir}/${prefix}.trimmed.derep.fasta -sizeout
# get lenghth distributions
length_stats ${outdir}/${prefix}.trimmed.derep.fasta  ${outdir}/${prefix}
# get stats for all the steps
seqkit stats ${outdir}/*${prefix}*.fa* > ${outdir}/${prefix}.stats
# run a fastqc
fastqc ${outdir}/${prefix}.3trimmed.fastq -o ${outdir}
# run the blast
#blast ${outdir}/${prefix}.trimmed.derep.fasta ${outdir}/${prefix} ${7}