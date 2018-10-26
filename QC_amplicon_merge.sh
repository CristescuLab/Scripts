#!/usr/bin/env bash
# usage: QC_amplicon_merge.sh primer_fwd primer_rev R1_filename
# for now only tested with CO1

qc(){
# To run function qc <fwd primer> <reverse primer> <R1 file> <num threads>
ADAPTER_FWD=$1
Adapter1rc=`echo ${ADAPTER_FWD}| tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
ADAPTER_REV=$2
Adapter2rc=`echo ${ADAPTER_REV}| tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
data=`basename ${3%%_R1.fastq.gz}`
data_folder=`dirname $3`
#runtrimmo
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 \
${data_folder}/${data}_R1.fastq.gz ${data_folder}/${data}_R2.fastq.gz \
${data}_R1_1P.fastq.gz ${data}_R1_1U.fastq.gz ${data}_R2_2P.fastq.gz \
${data}_R2_2U.fastq.gz ILLUMINACLIP:${adapters}:3:30:6 SLIDINGWINDOW:10:30

filterbytile.sh in1=${data}_R1_1P.fastq.gz in2=${data}_R2_2P.fastq.gz \
out1=${data}.R1.fastq.gz out2=${data}.R2.fastq.gz

cutadapt -j 10 -g "${ADAPTER_FWD}" -G "${ADAPTER_REV}" -a "${Adapter2rc}" \
-A "${Adapter1rc}" -o ${data}.1.fastq.gz -p ${data}.2.fastq.gz -m 200 \
--match-read-wildcards -q 28 --trim-n -l 250 ${data}.R1.fastq.gz \
${data}.R2.fastq.gz > ${data}.log1

pear -f ${outdir}/${data}.1.fastq.gz -r ${outdir}/${data}.2.fastq.gz \
-o ${outdir}/${data} -q 28 -t 150 -n 200 -m 600 -j 8 -e 2 -p 0.01 -v 20

seqkit -j 10 subseq -r 1:360 ${data}.assembled.fastq | seqkit -j 10 rmdup \
-s -i -o ${data}.clean.fastq -d ${pref}_duplicated.fastq.gz \
-D ${pref}_duplicated.detail.txt
}

# primers for 18S
f="GGWACWGGWTGAACWGTWTAYCCYCC"
r="TAAACTTCAGGGTGACCAAAAAATCA"
qc ${f} ${r} ${1} ${2}