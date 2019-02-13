#!/usr/bin/env bash
# index the reference
directory=/media/jshleap/ExtraDrive2/Playground/Mutation_accomulation
bwa index ${directory}/D_pulex_ref_PA42_clean.fasta -p PA42
for i in *.1.fastq.gz; do
    prefix=${i%%.1.fastq.gz}
    bwa mem -aHMpP -t 28 ${directory}/PA42 ${prefix}.1.fastq.gz  ${prefix}.2.fastq.gz | \
    samtools view -Shu - | samtools sort - | samtools rmdup -s - - | \
    tee ${prefix}_sorted_dedup.bam | bamToBed > ${prefix}_sorted_dedup.bed
done