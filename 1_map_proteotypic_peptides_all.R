source("R/general.R")
source("R/ConDigest.R")
library(stringr)

# download latest swissprot from uniprot
get_swissprot("swissprot-new") # runs automatically and will download a file with todays date. Insert that file in the field below

swissprot                   <- read.csv("data/20180503_swissprot(20341).csv", colClasses = rep("character", 20)) # change name
swissprot$Peptide           <- unlist(lapply(swissprot$Peptide, quest_del)) # get rid of questionmarks from columns. treat as they are modified
swissprot$Propeptide        <- unlist(lapply(swissprot$Propeptide, quest_del))
swissprot$Signal.peptide    <- unlist(lapply(swissprot$Signal.peptide, quest_del))

# Map modifications from swissprot
res <- list()

for(i in 1:dim(swissprot)[1]){
  pep           <- data.frame(entry = swissprot$Entry[i], ConDigest(swissprot[i,3], enzyme = "trypsin.strict")) # digest protein
  prot_mod      <- get_modifications(i) # get modifications
  prot_pep      <- get_peptide(i)
  prot_alt      <- get_altseq(i)
  pep$prot_mod  <- "" # empty column
  pep$prot_pep  <- "" # empty column
  pep$prot_alt  <- "" # empty column
  pep$max       <- swissprot$Length[i]
  
  # map protein modifications
  if(!length(prot_mod) == 1){
    for (j in 1:dim(prot_mod)[1]){
      filter <- check_interval(pep$fstart, pep$fstop, prot_mod[j, "start"], prot_mod[j, "stop"])
      pep[filter,"prot_mod"] <- paste0(pep[filter,"prot_mod"], paste(prot_mod$modification[j], "[", paste0(unique(c(prot_mod$start[j], prot_mod$stop[j])), collapse = ","), "]", sep =""), sep =";")
    }
  }
  
  # map protein peptide modifications
  if(!length(prot_pep) == 1){
    for (j in 1:dim(prot_pep)[1]){
      filter <- check_interval(pep$fstart, pep$fstop, prot_pep[j, "start"], prot_pep[j, "stop"])
      pep[filter,"prot_pep"] <- paste0(pep[filter,"prot_pep"], paste(prot_pep$class[j], "[", paste0(unique(c(prot_pep$start[j], prot_pep$stop[j])), collapse = ","), "]", sep =""), sep =";")
    }
  }
  
  # map protein alternative sequences plus known natural variants
  if(!length(prot_alt) == 1){
    for (j in 1:dim(prot_alt)[1]){
      filter <- check_interval(pep$fstart, pep$fstop, prot_alt[j, "start"], prot_alt[j, "stop"])
      pep[filter,"prot_alt"] <- paste0(pep[filter,"prot_alt"], paste(prot_alt$class[j], "[", paste0(unique(c(prot_alt$start[j], prot_alt$stop[j])), collapse = ","), "]", sep =""), sep =";")
    }
  }
  
  # filter out temp separators
  pep$prot_mod <- gsub("^;", "", pep$prot_mod)
  pep$prot_pep <- gsub("^;", "", pep$prot_pep)
  pep$prot_alt <- gsub("^;", "", pep$prot_alt)
  
  pep$prot_mod <- gsub(";$", "", pep$prot_mod)
  pep$prot_pep <- gsub(";$", "", pep$prot_pep)
  pep$prot_alt <- gsub(";$", "", pep$prot_alt)
  
  pep$len <- pep$stop-pep$start + 1
  res[[i]] <- pep[which(pep$len > 5 & pep$len <= 35),] # Apply relevant peptide filters
}

# rbind template
temp <- rbindlist(res)
temp$peptide <- paste(temp$peptide)
temp$fpeptide <- paste(temp$fpeptide)
write.csv(temp, paste0("output/", today, "_swissprot_peptides(", dim(temp)[1],").csv"), row.names = F)

