library(xlsx)
library(stringr)
library(data.table)
source("R/ConDigest.R")

df <- read.table("data/20160708.xls", 
                 sep = "\t", 
                 header = T, 
                 colClass = (c("character", "factor", "character", "numeric", "character", "character", "character", "factor", "factor")))

colnames(df) <- c("prest", "multi", "sequence", "weight", "gene", "ensg", "uniprot", "ms")
df <- df[df$ms == "GO", ] # filter out pretsts that have passed MS, figure 1e


prot    <- df[, c("prest", "sequence")]
digest  <- list()
for (i in 1:dim(prot)[1]){
  digest[[i]] <- data.frame(entry = prot[i,1], ConDigest(prot[i,2], missed = 0, enzyme = "trypsin.strict"))[, c("entry", "peptide")]
}
dig             <- data.frame(rbindlist(digest))
dig$entry       <- paste(dig$entry)
dig$peptide     <- paste(dig$peptide)

dig$nchar       <- unlist(lapply(paste(dig$peptide), nchar))
dig             <- data.table(dig[which(dig$nchar > 5 & dig$nchar < 35), ]) # get rid of short and longer peptides
prest_pep       <- dig[, coll(entry), by = peptide]

write.csv(prest_pep, paste0("output/20160608_prest(41119)_peptides(", length(unique(prest_pep$peptide)), ").csv"), row.names = F)
