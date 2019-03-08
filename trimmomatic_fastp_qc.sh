#!/usr/bin/env bash
# This is a test of Trimmomatic plus fastp to pseudoautomatically QC reads
# Arguments:
# 1) R1 file
# 2) CPUs to use
# 3) Required lenght after trimming
# 4) Trimming front
# 5) Minimum overlap
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
data=`basename ${1%%_R1.fastq.gz}`
data_folder=`dirname $1`
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 -threads $2 \
${data_folder}/${data}_R1.fastq.gz ${data_folder}/${data}_R2.fastq.gz \
${data}_R1_1P.fastq.gz ${data}_R1_1U.fastq.gz ${data}_R2_2P.fastq.gz \
${data}_R2_2U.fastq.gz ILLUMINACLIP:${adapters}:3:30:6 MINLEN:20


fastp -i ${data}_R1_1P.fastq.gz -I ${data}_R2_2P.fastq.gz -o ${data}.1.fastq.gz \
-O ${data}.2.fastq.gz --detect_adapter_for_pe --qualified_quality_phred 25 \
--length_required $3 --low_complexity_filter --overrepresentation_analysis \
--thread $2 --trim_front1 $4 --correction --overlap_len_require $5 -h ${data}

fastqc ${data}.*.fastq.gz

#pear -f /${data}.1.fastq.gz -r ${data}.2.fastq.gz -o ${data} -q 20 -t 150  \
#--min-overlap $5 --min-assembly-length $6 --max-assembly-length $7 -j $2 \
#--score-method 2 -p 0.01