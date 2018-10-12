# Title     : nucleotide diversity
# Objective : Compute nucleotide diversity in input alignments
# TAKEN AND MODIFIED FROM https://github.com/sophiadavid1/nucleotide_diversity_analysis
# Created by: jshleap
# Created on: 10/12/18

require(pegas)
require(ape)

file.names <- dir(pattern=".fasta")
args = commandArgs(trailingOnly=TRUE)

for(i in 1:length(file.names)){
	gene <- read.dna(file.names[i], format="fasta", as.matrix=FALSE)
	pi <- nuc.div(gene)
	write.table(file.names[i], file=args[1], append=TRUE)
	write.table(pi, file=args[1], append=TRUE)

}