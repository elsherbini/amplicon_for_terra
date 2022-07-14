library(dada2)
library(Biostrings)

seqtab.nochim <- readRDS(snakemake@input[["seqtab"]])

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

species_assignment <- assignSpecies(
                            dna,
                            snakemake@input[["silva_species_db"]],
                            allowMultiple = TRUE,
                            tryRC = TRUE,
                            n = 1000,
                            verbose = FALSE)


saveRDS(species_assignment, snakemake@output[["taxonomy_silva_species"]])

