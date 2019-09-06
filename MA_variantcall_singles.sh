#!/usr/bin/env bash
# exit when any command fails
set -e

#define reference genome path
#ref=/media/jshleap/ExtraDrive2/Playground/Mutation_accomulation/D_pulex_ref_PA42_clean.masked.fasta

# Usage:
# MA_variantcall_singles.sh path2reference input_prefix cpus
ref=${1}
java='java -jar -Xmx10g'
java='java -jar -Xmx10g'
GATK=~/Programs/gatk-4.1.0.0/GenomeAnalysisTK.jar
picard=~/Programs/picard/picard.jar
realSFS=~/Programs/angsd/misc/realSFS
density_plot(){
python3 - << EOF
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
df = pd.read_csv('$1', sep='\t')
fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for i, p in enumerate(['QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum']):
    ax = fig.add_subplot(2, 3, i+1)
    ax.title.set_text(p)
    d = df[p]
    m = float(d.mode().iloc[0])
    print(m)
    #d.plot.kde(bw_method=0.1, ax=ax)
    d.hist(ax=ax, density=True, bins=100)
    ax.axvline(x=m, lw=1, ls=':', color='b')
plt.tight_layout()
plt.savefig('$2_densities.pdf')
EOF
}


select_variants(){
${java} ${GATK} VariantsToTable -R ${2} -V ${1}_annotated.vcf \
-F ExcessHet -F AC -F AF -F AN -F DP -F FS -F MBQ -F QD -F SOR -F MQ \
-F MQRankSum -F ReadPosRankSum -O ${1}_annotations.table
modes=($(density_plot  ${1}_annotations.table test))
modes[1]=`bc -l <<< "scale=2;${modes[1]} +1"`
modes[2]=`bc -l <<< "scale=2;${modes[2]} +1"`
mq1=`bc -l <<< "scale=2;${modes[3]} + 0.1"`
mq2=`bc -l <<< "scale=2;${modes[3]} - 0.1"`
mqrs1=`bc -l <<< "scale=2;${modes[4]} - 0.5"`
mqrs2=`bc -l <<< "scale=2;${modes[4]} + 0.5"`
rprs1=`bc -l <<< "scale=2;${modes[5]} + 1.5"`
rprs2=`bc -l <<< "scale=2;${modes[5]} - 1.5"`
if [[ ! -f ${1}_SNPS.vcf ]]
then
    #Retain only biallelic and filter
    ${java} ${GATK} SelectVariants -V ${1}_raw_genos.vcf -O ${1}_SNPS.vcf \
    --restrict-alleles-to BIALLELIC -select-type SNP \
    -select "AF < 0.99" -select "SOR < ${modes[2]}" \
    -select "ReadPosRankSum < ${rprs1}" -select "ReadPosRankSum > ${rprs2}" \
    -select "MQRankSum > ${mqrs1}" -select "MQRankSum < ${mqrs2}" \
    -select "MQ > ${mq2}" -select "MQ < ${mq1}" \
    -select "FS < ${modes[1]}" -select "QD > ${modes[0]}"
fi

if [[ ! -f ${1}_SNP_db ]]
then
    # last selection step
    z --vcf ${1}_SNPS.vcf --max-missing 1 --maf 0.05 --minDP 20.0 \
    --min-meanDP 29.0 --max-meanDP 31 --minQ 50.0 --recode --recode-INFO-all \
    --min-alleles 2 --max-alleles 2 --hwe 0.001 --out ${1}_SNP_db
fi
}

filter_bam(){
# 1. ref
# 2. bam file

if [[ ! -f ${2}.bai ]]; then
    samtools index -@ 28 $2
fi
lab=${2%%_markdup.bam}

if [[ ! -f ${lab}_raw_variants.vcf ]]
then
    ${java} ${GATK} HaplotypeCaller -R $1 -I $2 -ERC GVCF \
    -O ${lab}_raw_variants.vcf -A StrandOddsRatio -A ReadOrientationArtifact \
    -A AlleleFraction -A BaseQuality -A BaseQualityRankSumTest -A Coverage \
    -A DepthPerAlleleBySample -A DepthPerSampleHC -A LikelihoodRankSumTest \
    -A MappingQuality -A MappingQualityRankSumTest -A ExcessHet \
    -A QualByDepth  -A StrandBiasBySample  -A FisherStrand \
    --standard-min-confidence-threshold-for-calling 30
fi

if [[ ! -f ${lab}_raw_genos.vcf ]]
then
    ${java} ${GATK} GenotypeGVCFs -R ${ref} -V ${lab}_raw_variants.vcf \
    -O ${lab}_raw_genos.vcf
fi

if [[ ! -f ${lab}_annotated.vcf ]]
then
    ${java} ${GATK} VariantAnnotator -R ${ref} -I ${lab}_markdup.bam \
   -V ${lab}_raw_genos.vcf -O ${lab}_annotated.vcf -A StrandOddsRatio \
   -A ReadOrientationArtifact -A AlleleFraction -A BaseQuality \
   -A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample \
   -A DepthPerSampleHC -A LikelihoodRankSumTest -A MappingQuality \
   -A MappingQualityRankSumTest -A ExcessHet -A QualByDepth \
   -A StrandBiasBySample -A FisherStrand
fi
# select_variants ${lab}
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
tag=${3%%_markdup.bam}
for i in `seq $1`
do
   echo "Running iteration $i in bootstrap"
   select_variants ${tag} ${ref}
   ${java} ${GATK} IndexFeatureFile -F ${tag}_SNP_db.recode.vcf
   ${java} ${GATK} BaseRecalibrator -R ${ref} -I ${tag}_markdup.bam \
   -O ${tag}_recal_data.table --known-sites ${tag}_SNP_db.recode.vcf
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
angsd -bam ${tag}_markdup.bam -doSaf 1 -anc ${ref} -GL 1 -P ${3} -out ${tag}.out
#angsd -bam <(ls ${tag}_markdup.bam) -doSaf 1 -anc ${ref} -GL 1 -P ${3} \
#-out ${tag}.out
${realSFS} out.saf.idx >${tag}_est.ml

ehe=`python - << EOF
import numpy
a = numpy.loadtxt('${tag}_est.ml')
print(a[1]/a.sum())
EOF`

