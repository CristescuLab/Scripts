#!/usr/bin/env bash
# exit when any command fails
set -e

#define reference genome path
ref=/media/jshleap/ExtraDrive2/Playground/Mutation_accomulation/D_pulex_ref_PA42_clean.fasta.masked


filter_bam(){
# 1. ref
# 2. bam file
java='java -jar -Xmx10g'
GATK=~/Programs/gatk-4.1.0.0/GenomeAnalysisTK.jar
picard=~/Programs/picard/picard.jar

if [[ ! -f ${2}.bai ]]; then
    samtools index -@ 28 $2
fi

${java} ${GATK} HaplotypeCaller -R $1 -I $2 -ERC GVCF -O raw_variants.vcf
--standard-min-confidence-threshold-for-calling 30 \
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
for i in seq $1; do
   ${java} ${GATK} IndexFeatureFile -F SNP_db.recode.vcf
   ${java} ${GATK} BaseRecalibrator -R ${ref} -I ${tag}_markdup.bam \
   -O ${tag}_recal_data.table --known-sites SNP_db.vcf
done

}

# index reference
bwa index -a bwtsw ${ref}
# create dictionary
picard=~/Programs/picard/picard.jar
${java} ${picard} CreateSequenceDictionary R=$2 O=${2%%fasta}.dict
for i in *.1.fastq.gz; do
    # get name of pair
    tag=${i%%.1.fastq.gz}
    # map files to reference and ignore unmapped reads (this is to avoid possible contaminant sequences)
    bwa mem -t 28 -aM -R "@RG\tID:${tag}\tSM:${tag}\tPL:Illumina" ${ref} \
    ${tag}.1.fastq.gz ${tag}.2.fastq.gz | samtools view -b -h -F 4 - > ${tag}_mapped.bam
    # Sort mapped files
    java -jar ~/Programs/picard/picard.jar SortSam I=${tag}_mapped.bam \
    O=${tag}_sorted.bam SORT_ORDER=coordinate
    # Mark duplicates and remove them
    java -jar ~/Programs/picard/picard.jar MarkDuplicates I=${tag}_sorted.bam \
    O=${tag}_markdup.bam M=${tag}_dup_metrics.txt \
    REMOVE_SEQUENCING_DUPLICATES=true ASSUME_SORTED=true
    # Get some summary stats
    java -jar ~/Programs/picard/picard.jar CollectAlignmentSummaryMetrics \
    R=${ref} I=${tag}_markdup.bam OUTPUT=C001-88.stats
    # get depth at each position
    samtools depth ${tag}_markdup.bam > ${tag}.depth
    samtools index ${tag}_markdup.bam
    # boostrap


done





# compute expected heterocigocity
angsd -bam <(ls *_markdup.bam) -doSaf 1 -anc ${ref} -GL 1 -P 24 -out out
~/Programs/angsd/misc/realSFS out.saf.idx >est.ml

ehe=`python - << EOF
import numpy
a = numpy.loadtxt('est.ml')
print(a[1]/a.sum())
EOF`


# Merge files
finp=`ls -1 *_markdup.bam| paste -sd "," -| sed 's/,/ I=/g'
inp=`echo " I=${finp}"`
java -jar ~/Programs/picard/picard.jar MergeSamFiles ${inp} O=output_merged_files.bam USE_THREADING=true

# Add the ancestor orphaned read group for accomulate
samtools addreplacerg -r "ID:ancestor" -r "SM:ancestor" -m orphan_only  output_merged_files.bam > CAN.bam

# Sort the resulting bam file
java -jar ~/Programs/picard/picard.jar SortSam I=CAN.bam o=CAN_sorted.bam SORT_ORDER=coordinate

#Run accumulate
samtools view -@ 28 -H CAN_sorted.bam| ~/Programs/accuMulate-tools/extract_samples.py ancestor - > params.ini
~/Programs/accuMulate-tools/GC_content.py ${ref} >> params.ini
echo "mu=2.30e-9
seq-error=0.001
phi-diploid=0.001
phi-haploid=0.001
ploidy-ancestor=2
ploidy-descendant=2" >> params.ini

~/Programs/accuMulate-tools/dictionary_converter.py ${ref} > reference_genome.dict
mkdir -p tmp
bedtools makewindows -g reference_genome.dict -w 100000 | split -l 1 - tmp/
parallel -j 28 accuMUlate --theta ${ehe} -c params.ini -b CAN_sorted.bam \
-r ${ref} -i {} ::: tmp/* > unsorted_out.txt
