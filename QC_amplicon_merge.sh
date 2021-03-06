#!/usr/bin/env bash
# usage: QC_amplicon_merge.sh primer_fwd primer_rev R1_filename outputdir
# for now only tested with CO1
# If used on computecanada:
# module load trimmomatic fastqc intel/2016.4 bbmap/37.36 pear/0.9.10

qc(){
# create output folder if required
mkdir -p $4
# To run function qc <fwd primer> <reverse primer> <R1 file> <num threads>
ADAPTER_FWD=$1
Adapter1rc=`echo ${ADAPTER_FWD}| tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
ADAPTER_REV=$2
Adapter2rc=`echo ${ADAPTER_REV}| tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
data=`basename ${3%%_R1.fastq.gz}`
data_folder=`dirname $3`

# Run fastqc on the raw reads
fastqc ${data_folder}/${data}_R[12].fastq.gz -o ${4}
#runtrimmo
adapters=${EBROOTTRIMMOMATIC}/adapters/NexteraPE-PE.fa
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.36.jar PE -phred33 \
${data_folder}/${data}_R1.fastq.gz ${data_folder}/${data}_R2.fastq.gz \
${4}/${data}_R1_1P.fastq.gz ${4}/${data}_R1_1U.fastq.gz \
${4}/${data}_R2_2P.fastq.gz ${4}/${data}_R2_2U.fastq.gz \
ILLUMINACLIP:${adapters}:3:30:6 SLIDINGWINDOW:10:20

# Run fastqc on the trimmomatic output
fastqc ${4}/${data}_R[12]_[12][UP].fastq.gz -o ${4}

# run filter tile
filterbytile.sh in1=${4}/${data}_R1_1P.fastq.gz in2=${4}/${data}_R2_2P.fastq.gz \
out1=${4}/${data}.R1.fastq.gz out2=${4}/${data}.R2.fastq.gz usejni=t

# run fastqc on the filter tile
fastqc ${4}/${data}.R*.fastq.gz -o ${4}

# remove primers with cutadappt
cutadapt -j 10 -g "${ADAPTER_FWD}" -G "${ADAPTER_REV}" -a "${Adapter2rc}" \
-A "${Adapter1rc}" -o ${4}/${data}.1.fastq.gz -p ${4}/${data}.2.fastq.gz \
-m 200 --match-read-wildcards -q 20 --trim-n -l 250  ${4}/${data}.R1.fastq.gz \
 ${4}/${data}.R2.fastq.gz > ${data}.log1

# run fastqc on the Cutadapt output
fastqc ${4}/${data}.[12].fastq.gz -o ${4}

pear -f ${4}/${data}.1.fastq.gz -r ${4}/${data}.2.fastq.gz \
-o ${4}/${data} -q 20 -t 150 -n 200 -m 600 -j 8 -e 2 -p 0.01 -v 20

# run fastqc on the pear  output
fastqc ${4}/${data}*assembled*.fastq -o ${4}

# Subset and deduplicate
seqkit -j 10 subseq -r 1:360 ${4}/${data}.assembled.fastq | seqkit -j 10 rmdup \
-s -i -o ${4}/${data}.clean.fastq -d ${4}/${data}_duplicated.fastq.gz \
-D ${4}/${data}_duplicated.detail.txt

seqkit -j 10 fq2fa ${4}/${data}.clean.fastq -o ${4}/${data}.clean.fasta

# deduplicate unasembled forwardf
seqkit -j 10 rmdup -s -i -o ${4}/${data}.unassembled.forward.clean.fastq \
-d ${data}.unassembled.forward.duplicated.fastq.gz \
-D ${data}.unassembled.forward.duplicated.detail.txt \
${4}/${data}.unassembled.forward.fastq

seqkit -j 10 fq2fa ${4}/${data}.unassembled.forward.clean.fastq \
-o ${4}/${data}.unassembled.forward.clean.fasta

# deduplicate unasembled reverse
seqkit -j 10 rmdup -s -i -o ${4}/${data}.unassembled.reverse.clean.fastq \
-d ${4}/${data}.unassembled.reverse.duplicated.fastq.gz \
-D ${4}/${data}.unassembled.reverse.duplicated.detail.txt \
${4}/${data}.unassembled.reverse.fastq

seqkit -j 10 fq2fa ${4}/${data}.unassembled.reverse.clean.fastq \
-o ${4}/${data}.unassembled.reverse.clean.fasta

# run fastqc on the last output
fastqc ${4}/${data}.clean.fastq -o ${4}
}

# primers for 18S
#f="AGGGCAAKYCTGGTGCCAGC"
#r="GRCGGTATCTRATCGYCTT"
qc ${1} ${2} ${3} ${4}