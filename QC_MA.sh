#!/usr/bin/env bash

adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
data=`basename ${1%%_R1.fastq.gz}`
data_folder=`dirname $1`
#runtrimmo
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 -threads $2 \
${data_folder}/${data}_R1.fastq.gz ${data_folder}/${data}_R2.fastq.gz \
${data}_R1_1P.fastq.gz ${data}_R1_1U.fastq.gz ${data}_R2_2P.fastq.gz \
${data}_R2_2U.fastq.gz ILLUMINACLIP:${adapters}:3:30:6 SLIDINGWINDOW:10:30 \
MINLEN:20

filterbytile.sh in1=${data}_R1_1P.fastq.gz in2=${data}_R2_2P.fastq.gz \
out1=${data}.R1.fastq.gz out2=${data}.R2.fastq.gz

cutadapt -j $2 -a "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-A "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-G "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -o ${data}.1.fastq.gz \
-p ${data}.2.fastq.gz --match-read-wildcards -q 28,28 --trim-n -m 20 \
 ${data}.R1.fastq.gz ${data}.R2.fastq.gz > ${data}.log1

fastqc ${data}.[12].fastq.gz