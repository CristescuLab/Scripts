# Title     : nucleotide diversity
# Objective : Compute nucleotide diversity in input alignments
# TAKEN AND MODIFIED FROM https://github.com/sophiadavid1/nucleotide_diversity_analysis
# Created by: jshleap
# Created on: 10/12/18
# In cedar need to load:
# module load gcc/5.4.0 r-bundle-bioconductor/3.4

require(pegas)
require(ape)

file.names <- dir(pattern=".fas$")
args = commandArgs(trailingOnly=TRUE)

df <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("File", "pi")
colnames(df) <- x
j=1
for(i in 1:length(file.names)){
	gene <- read.dna(file.names[i], format="fasta", as.matrix=FALSE)
    if (length(gene) > 2){
	pi <- nuc.div(gene)
	df[j,"File"] = file.names[i]
    df[j,"pi"] = pi
    j = j+1}
}
write.table(df, args[1], sep='\t')