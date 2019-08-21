# Title     : TODO
# Objective : TODO
# Created by: jshleap
# Created on: 2019-08-20
suppressPackageStartupMessages(library(dada2));
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(phylotools))

option_list <- list(
make_option(c("-p", "--pattern"), type="character", default=NULL,
help="Pattern where to the dechimerixed (nochimaeras.rda) files"),
make_option(c("-o", "--outprefix"), type="character", default='all_samples',
help="Prefix to ouputs")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

pattern <- opt$pattern
prefix <- opt$outprefix
rdas <- Sys.glob(pattern)
rda_list <- list()
count <- 1
for (rd in rdas) {
    load(rd)
    rda_list[[count]] <- seqtab.nochim
    count <- count + 1
    rm(seqtab.nochim)
}

fn <- paste0(prefix, '.fasta')
merged_table <- mergeSequenceTables(tables=rda_list)
all <- getUniques(merged_table)
uniquesToFasta(all, fn)
seqnames <- get.fasta.name(paste0(prefix, '.fasta'), clean_name = FALSE)
otutab <- merged_table
colnames(otutab) <- seqnames
write.table(otutab, file=paste0(prefix, '_ASV_table.tsv'), sep='\t')
dir.create('ASVs_MERGED')
for(sample in row.names(merged_table)){
    print(sample)
    st <- merged_table[sample,]
    boolean <- st > 0
    st <- st[boolean]
    if (length(st) != 0){
        seqs <- DNAStringSet(getSequences(st))
        names(seqs) <- seqnames[boolean]
        outfn <- file.path("./ASVs_MERGED", paste0(sample,'_ASV.fasta'))
        writeXStringSet(seqs, outfn)
    }
}