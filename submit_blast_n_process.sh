#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -A def-mcristes
#SBATCH --job-name=467_2S

module purge
module load fastqc
module load nixpkgs/16.09
module load gcc/5.4.0
module load trimmomatic
module load fastx-toolkit/0.0.14
module load nixpkgs/16.09  gcc/5.4.0 blast+
module load scipy-stack/2018b

## This is the basic submission for blast and postprocess it.
cd $SLURM_SUBMIT_DIR

#have a nice day! :)
# The job command(s):
out=`cut -d'.' -f5 ${file_path}`
path="${file_path}/${out}"

blastn -db ${db_path} -query ${file_path}/*.trimmed.derep.fasta \
-evalue 0.0001 -perc_identity 97 -max_target_seqs 40 \
-outfmt "6 qseqid sseqid pident evalue qcovs stitle" -out ${path}.hits -num_threads 8
python blast.processing.py ${path}.hits ${path}_98_98 -p 98 -q 98
python blast.processing.py ${path}.hits ${path}_98_90 -p 98 -q 90
