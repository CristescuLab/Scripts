#!/usr/bin/env bash
#----------------------------------------------------------------------------------------------------------
# Code developed by Gillian Martin and Frederic Chain, 2017.
# Modified by Joanne Littlefair, and Jose Sergio Hleap, 2018.
#----------------------------------------------------------------------------------------------------------
# Usage:
# bash QC_pipeline.sh <path2fileR1> <Fprimer> <Rprimer>
# requires 8 cpus-per-task

# Set working directory where you submitted the job
cd $SLURM_SUBMIT_DIR

#Set variables and directories
## Input arguments
data_folder=`dirname $1`
data=`basename $1`
data=${data%%_R1*}
outdir=output_${data}
Fprimer=$2
Rprimer=$3
## Reverse complement sequences in bash:
Adapter1rc=`echo ${Fprimer} | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
Adapter2rc=`echo ${Rprimer} | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
echo "The forward primer sequence used is $Fprimer, and reverse primer $Rprimer"

## Paths to software
flash=~/projects/def-mcristes/Software/FLASH-1.2.11/flash
## Make output directory if it doesn't already exist, change directory name for each region
if [[ ! -d ${outdir} ]]; then
  mkdir -p ${outdir}
fi

#Sort sequences into reads with your primers allowing up to 2 missmatches
m=2 # set more if you need to (not recommended)
for i in R1 R2;do
seqkit -j 8 grep --pattern-file <(zcat ${data}_${i}.fastq.gz | \
seqkit -j 8 locate -i -d -p -m ${m} "${Fprimer},${Rprimer}"| \
awk '{if(NR>1)print $1}') ${data}_${i}.fastq.gz| \
seqkit -j 8 rmdup -n -i -o ${outdir}/${data}_sorted_${i}.fastq.gz
done

# Trim illumina adapters with trimmomatic

adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa # change if different adapter file
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 \
${outdir}/${data}_sorted_R1.fastq.gz ${outdir}/${data}_sorted_R2.fastq.gz \
${outdir}/${data}_R1_1P.fastq.gz ${outdir}/${data}_R1_1U.fastq.gz \
${outdir}/${data}_R2_2P.fastq.gz ${outdir}/${data}_R2_2U.fastq.gz \
ILLUMINACLIP:${adapters}:3:30:6 SLIDINGWINDOW:10:20

# Merge the reads with a pvalue of 0.01 for the test (-p 0.01; valid values:  0.0001, 0.001, 0.01, 0.05 and 1.0)
# Minimum overlap of 10 bp (-v 10), minimum lenght of assemby of 100 (-n 100),
# Trimming reads with less than 20 in consecutive quality (-q 20)
# with at least 100 bp per read after quality trimming (-t 100), running on 8 cpus (-j 8)
# and weight the assembly with the quality score
pear -f ${outdir}/${data}_R1_1P.fastq.gz -r ${outdir}/${data}_R2_2P.fastq.gz \
-o ${outdir}/${data} -q 20 -t 100 -n 100 -j 8 -e 2 -p 0.01 -v 10


# Remove trailing adapter / primer. This also discards reads less than 100bp into a file called ${data}.short.fastq
# On merged reads we use revcomp reverse primer and usual forward primer
cutadapt -a $Adapter2rc -m 100 --too-short-output ${outdir}/${data}.short.fastq \
-o ${outdir}/${data}.3trimmed.fastq ${outdir}/${data}.assembled.fastq
cutadapt -g ^$Fprimer -o ${outdir}/${data}.5trimmed.fastq ${outdir}/${data}.3trimmed.fastq

### Rerun fastqc on this filtered file
fastqc ${outdir}/${data}.5trimmed.fastq -o ${outdir}

## Dereplicate to reduce complexity within a single library (makes it easier to create region wide ZOTUs).
## This creates a report of what is being removed
seqkit -j 8 rmdup -s -i -m -o ${outdir}/${data}.trimmed.derep.fasta \
-d ${outdir}/${data}_duplicated.fa.gz -D ${outdir}/${data}_duplicated.detail.txt \
${outdir}/${data}.5trimmed.fastq
