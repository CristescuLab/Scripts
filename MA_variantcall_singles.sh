#!/usr/bin/env bash
# exit when any command fails
set -e

#define reference genome path
#ref=/media/jshleap/ExtraDrive2/Playground/Mutation_accomulation/D_pulex_ref_PA42_clean.masked.fasta
ref=${1}
java='java -jar -Xmx10g'

filter_bam(){
# 1. ref
# 2. bam file
java='java -jar -Xmx10g'
GATK=~/Programs/gatk-4.1.0.0/GenomeAnalysisTK.jar
picard=~/Programs/picard/picard.jar

if [[ ! -f ${2}.bai ]]; then
    samtools index -@ 28 $2
fi
lab=${2%%_markdup.bam}
${java} ${GATK} HaplotypeCaller -R $1 -I $2 -ERC GVCF \
-O ${lab}_raw_variants.vcf --standard-min-confidence-threshold-for-calling 30 \
-A AlleleFraction -A BaseQuality -A BaseQualityRankSumTest -A Coverage \
-A DepthPerAlleleBySample -A DepthPerSampleHC -A LikelihoodRankSumTest \
-A MappingQuality -A MappingQualityRankSumTest -A ExcessHet -A FisherStrand \
-A QualByDepth -A ReadOrientationArtifact -A StrandBiasBySample \
-A StrandOddsRatio

${java} ${GATK} GenotypeGVCFs -R ${ref} -V raw_variants.vcf -O raw_genos.vcf


#Retain only biallelic and fileter
${java} ${GATK} SelectVariants -V raw_genos.vcf -O SNPS.vcf \
--restrict-alleles-to BIALLELIC -select-type SNP \
-select "AF < 0.99" -select "SOR < 1.0" -select "ReadPosRankSum < 1.5" \
-select "ReadPosRankSum > 0.0" -select "MQRankSum > -0.5" \
-select "MQRankSum < 1.5" -select "MQ > 59.9" -select "MQ < 60.1" \
-select "FS < 6.0" -select "QD > 10.0"

# last selection step
vcftools --vcf SNPS.vcf --max-missing 1 --maf 0.05 --minDP 20.0 \
--min-meanDP 29.0 --max-meanDP 31 --minQ 50.0 --recode --recode-INFO-all \
--min-alleles 2 --max-alleles 2 --hwe 0.001 --out SNP_db
}

bootsrap_vcf(){
# run un-calibrated data to vcf, filter and repeat
# 1) number of repeats
# 2) reference
# 3) un-calibrated data
# first call
java='java -jar -Xmx10g'
GATK=~/Programs/gatk-4.1.0.0/GenomeAnalysisTK.jar
filter_bam $2 $3
for i in `seq $1`
do
   echo "Running iteration $i in bootstrap"
   ${java} ${GATK} IndexFeatureFile -F SNP_db.recode.vcf
   ${java} ${GATK} BaseRecalibrator -R ${ref} -I ${tag}_markdup.bam \
   -O ${tag}_recal_data.table --known-sites SNP_db.vcf
   echo ${i} > bootdone
done
}

# get name of pair
tag=${2} # prefix is a command line argument
# map files to reference and ignore unmapped reads (this is to avoid possible contaminant sequences)
if [[ ! -f  ${tag}_mapped.bam ]]
then
    bwa mem -t 28 -aM -R "@RG\tID:${tag}\tSM:${tag}\tPL:Illumina" ${ref} \
    ${tag}.1.fastq.gz ${tag}.2.fastq.gz | samtools view -b -h \
    -F 4 - > ${tag}_mapped.bam
fi
# Sort mapped files
if [[ ! -f ${tag}_sorted.bam ]]
then
    java -jar ~/Programs/picard/picard.jar SortSam I=${tag}_mapped.bam \
    O=${tag}_sorted.bam SORT_ORDER=coordinate
fi
# Mark duplicates and remove them
if [[ ! -f ${tag}_markdup.bam ]]
then
    java -jar ~/Programs/picard/picard.jar MarkDuplicates \
    I=${tag}_sorted.bam O=${tag}_markdup.bam M=${tag}_dup_metrics.txt \
    REMOVE_SEQUENCING_DUPLICATES=true ASSUME_SORTED=true
fi
# Get some summary stats
if [[ ! -f ${tag}.stats ]]
then
    java -jar ~/Programs/picard/picard.jar CollectAlignmentSummaryMetrics \
    R=${ref} I=${tag}_markdup.bam OUTPUT=${tag}.stats
fi
# get depth at each position
samtools depth ${tag}_markdup.bam > ${tag}.depth
samtools index ${tag}_markdup.bam
# boostrap
if [[ ! -f bootdone ]]
then
    bootsrap_vcf 100 ${ref} ${tag}_markdup.bam
else
    n=`cat bootdone`
    n=$(( 100 - n ))
    bootsrap_vcf ${n} ${ref} ${tag}_markdup.bam
fi

# compute expected heterocigocity
angsd -bam <(ls ${tag}_markdup.bam) -doSaf 1 -anc ${ref} -GL 1 -P 24 \
-out ${tax}.out
~/Programs/angsd/misc/realSFS out.saf.idx >${tag}_est.ml

ehe=`python - << EOF
import numpy
a = numpy.loadtxt('${tag}_est.ml')
print(a[1]/a.sum())
EOF`

