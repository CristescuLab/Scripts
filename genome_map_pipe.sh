#!/usr/bin/env bash
# This script will use the principles of Heidi Lischer for a reference-guided
# genome assembly, but modifying a couple of steps
# Jose Sergio Hleap 2018
# Assumes quality trimmed files



QC(){
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
data=`basename ${1%%_R1.fastq.gz}`
data_folder=`dirname $1`
#runtrimmo
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 \
${data_folder}/${data}_R1.fastq.gz ${data_folder}/${data}_R2.fastq.gz \
${data}_R1_1P.fastq.gz ${data}_R1_1U.fastq.gz ${data}_R2_2P.fastq.gz \
${data}_R2_2U.fastq.gz ILLUMINACLIP:${adapters}:3:30:6 SLIDINGWINDOW:10:30 \
MINLEN:40
# Deal with some of the tile issue
filterbytile.sh in1=${data}_R1_1P.fastq.gz in2=${data}_R2_2P.fastq.gz \
out1=${data}.R1.fastq.gz out2=${data}.R2.fastq.gz
# Run filtertile in the unpaired
filterbytile.sh in=${data}_R1_1U.fastq.gz out=${data}_R1_unpaired.fastq.gz
filterbytile.sh in=${data}_R2_2U.fastq.gz out=${data}_R2_unpaired.fastq.gz
# Remove adapters in pair Just incase something missing in trimmo
cutadapt -a "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-A "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-G "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -o ${data}.1.fastq.gz \
-p ${data}.2.fastq.gz --match-read-wildcards -q 28,28 --trim-n -m 40 \
 ${data}.R1.fastq.gz ${data}.R2.fastq.gz > ${data}.log1
# Remove adapters in unpaired Just in case something missing in trimmo
cutadapt -a "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" -a "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" \
-g "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -o ${data}.1U.fastq.gz \
--match-read-wildcards -q 28,28 --trim-n -m 40 \
${data}_R1_unpaired.fastq.gz > ${data}.log2
cutadapt -a "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
-g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" -a "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" \
-g "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" -o ${data}.2U.fastq.gz \
--match-read-wildcards -q 28,28 --trim-n -m 40 \
${data}_R2_unpaired.fastq.gz > ${data}.log3
cat ${data}.*U.fastq.gz > ${data}.unpaired.fastq.gz
# Quality
fastqc ${data}*.fastq.gz
}

define_blocks(){
sample=$1 # sample name
line=`cut -d'-' -f1 ${sample}`
ref=PA42  # give it a reference database name
query=$3 # query prefix
reffas=$4 # reference fasta file name
#remove scaffolds shorter than 10 kb
cat ${reffas}| seqkit seq -V 0 -m 10000 > ${reffas%%.fasta}_clean.fasta
reffas=${reffas%%.fasta}_clean.fasta
#BWA, create db first
if [ ! -f PA42.bwt ];then
    bwa index -p ${ref} -a bwtsw ${reffas}
fi

if [ ! -f ${sample}_all.sorted.bam ]; then
    # Map paired
    bwa mem -aHM -R "@RG\tID:${line}\tSM:${line}" ${ref} ${query}.1.fastq.gz \
    ${query}.2.fastq.gz > ${sample}.bwa.sam
    # Map unpaired
    bwa mem -aHM -R "@RG\tID:${line}\tSM:${line}" ${ref} \
    ${query}.unpaired.fastq.gz >> ${sample}.bwa.sam
    cat ${sample}.bwa.sam | samtools view -bS - | samtools sort - -T ${sample} \
    -o ${sample}_all.sorted.bam
    rm ${sample}.bwa.sam
    samtools index ${sample}_all.sorted.bam
fi

#filter mapped reads
samtools view -b F 4 ${sample}_all.sorted.bam > ${sample}.sorted.bam
samtools index ${sample}.sorted.bam

#get unmapped reads
samtools view -b -f 4 ${sample}_all.sorted.bam > ${sample}_unmapped.sorted.bam
samtools view -b -f 9 ${sample}_unmapped.sorted.bam | bedtools bamtofastq \
-fq ${sample}_failPair.1.fastq -fq2 ${sample}_failPair.2.fastq
samtools view -b -F 8 -f 64 ${sample}_unmapped.sorted.bam | bamtools convert \
-format fastq -out ${sample}_failUnPairR1.fastq
samtools view -b -F 8 -f 128 ${sample}_unmapped.sorted.bam | bamtools convert \
-format fastq -out ${sample}_failUnPairR2.fastq
bamtools stats -in  _all.sorted.bam >> ${sample}_blocks.log

#filter for mapping quality >=10
samtools view -b -q 10 ${sample}.sorted.bam > ${sample}.sorted.filtered.bam
bamtools stats -in  >> ${sample}.sorted.filtered.bam >> ${sample}_blocks.log
echo "get blocks and superblocks..." >> ${sample}_blocks.log
covFile=${sample}_coverage.txt
}

run_anviomap(){
sample=$1 # sample name
ref=$2 # give it a reference database nam
query=$3 # query fasta file
reffas=$4 # reference fasta file name
datadir=$5 # where the original fastq file stored
output_dir=$6 #final output directory


# index it in a anvio-friendly manner
anvi-init-bam ${sample}.bwa-raw.bam -o ${sample}

}
