library(DECIPHER)
library(dada2)
library(Biostrings)

seqtab.nochim <- readRDS(snakemake@input[["seqtab"]])

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(snakemake@input[["decipher_db"]]) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand = "both", processors = snakemake@threads, verbose = FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)

saveRDS(taxid, snakemake@output[["taxonomy_decipher"]])
