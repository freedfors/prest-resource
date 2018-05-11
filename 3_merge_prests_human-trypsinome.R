# Denna som används för att producera pltos till publikationen 

library(data.table)
source("R/UniquePep.R")

#

prest_pep <- data.table(read.csv("output/20160608_prest(41119)_peptides(205822).csv")) # all prests and their peptides
#all_pep   <- read.csv("output/20160708_swissprot_peptides(627331).csv") 
all_pep   <- read.csv("output/20180503_swissprot_peptides(629341).csv")

# eliminate duplicated peptides from list
all_pep$peptide <- paste(all_pep$peptide) # as character
all_pep$uni     <- unlist(lapply(paste(all_pep$peptide), UniquePep)) # unique peptides based on swissprot in-silico digestion
uni_pep         <- all_pep[all_pep$uni == TRUE,]

combined  <- merge(prest_pep, uni_pep, by.y = "peptide", by.x = "peptide", all.y = T, all.x = F)
combined  <- combined[which(!is.na(combined$peptide)),] # get rid of NA orginating from prest_pep
combined  <- combined[which(!duplicated(combined))] # get rid of duplicated rows, including start stop
#saveRDS(combined, "data/20160608_peptides_all_prest_combined.RDa")
#combined <- readRDS("data/20160608_peptides_all_prest_combined.RDa")
#combined <- readRDS("data/20160608_peptides_all_prest_combined.RDa")

# barplot for evaluation
combined$before3 <- unlist(lapply(combined$fpeptide, function(x){substring(x, 1,1)}))
combined$before2 <- unlist(lapply(combined$fpeptide, function(x){substring(x, 2,2)}))
combined$before1 <- unlist(lapply(combined$fpeptide, function(x){substring(x, 3,3)}))
combined$after3 <- unlist(lapply(paste(combined$fpeptide), function(x){aa_after(x, shift = 0)}))
combined$after2 <- unlist(lapply(paste(combined$fpeptide), function(x){aa_after(x, shift = 1)}))
combined$after1 <- unlist(lapply(paste(combined$fpeptide), function(x){aa_after(x, shift = 2)}))

########################################
# ugly code but fastest solution atm ###
########################################

# start early
#fstart = -2
combined[which(combined$fstart == -2), c("before3", "before2", "before1")] <- NA
# fstart = -1
combined[which(combined$fstart == -1), c("before1")] <- combined[which(combined$fstart == -1), c("before3")]
combined[which(combined$fstart == -1), c("before3", "before2")] <- NA
# fstart
combined[which(combined$fstart == 0), c("before2", "before1")] <- combined[which(combined$fstart == 0), c("before3", "before2")]
combined[which(combined$fstart == 0), c("before3")] <- NA

# end of peptide sequence
combined$late <- combined$max - combined$fstop
# late = -1
combined[which(combined$late == -1), c("after2", "after1")] <- combined[which(combined$late == -1), c("after3", "after2")]
combined[which(combined$late == -1), c("after3")] <- NA
# late = -2
combined[which(combined$late == -2), c("after1")] <- combined[which(combined$late == -2), c("after3")]
combined[which(combined$late == -2), c("after2" , "after3")] <- NA
# late = -3
combined[which(combined$late == -3), c("after1", "after2" , "after3")] <- NA
####

combined$first        <- substring(combined$peptide, 1,1)
combined$second       <- substring(combined$peptide, 2,2)
combined$third        <- substring(combined$peptide, 3,3)
combined$last         <- unlist(lapply(combined$peptide, function(x){substring(x, nchar(x), nchar(x))}))
combined$last_second  <- unlist(lapply(combined$peptide, function(x){substring(x, (nchar(x)-1), (nchar(x)-1))}))
combined$last_third   <- unlist(lapply(combined$peptide, function(x){substring(x, (nchar(x)-2), (nchar(x)-2))}))
combined$met          <- unlist(lapply(combined$peptide, function(x){str_count(x, "M")}))

# function to calculate amino acid after
aa_after <- function(x, shift = 0){
  substring(x, (nchar(x)-shift), (nchar(x)- shift))
}


#####

#saveRDS(combined, "data/20180503_peptides_all_prest_combined_flanking.RDa")
combined <- readRDS("data/20180503_peptides_all_prest_combined_flanking.RDa")
#######

# visualize trypsinome
combined <- data.frame(combined) # create data frame

if (dim(combined[which(duplicated(combined$peptide)),])[1] > 0){ # control for duplicated peptides, must be
  warning("Duplicated peptides")
}
combined <- combined[!duplicated(combined$peptide),]
if(length(which(combined$uni == F) > 0)){ # control that only gene specific peptides are included
  warning("Non-unique peptides")
}

combined  <-combined[which(combined$last %in% c("R", "K") & combined$before1 %in% c("R", "K")),] #only tryptic peptide
plot1     <- table(combined$len) # only tryptic peptides
sub1      <- combined[which(combined$prot_mod == "" & combined$prot_pep == ""),] # only non modified residues
plot2     <- table(sub1$len) # peptides wo glycogroups +/- peptides away
sub1a     <- sub1[which(sub1$met == 0),]
plot2a    <- table(sub1a$len)
#sub2 <- sub1[which(sub1$prot_alt == ""),]
#plot3 <- table(sub2$len)
# Filter out peptides with DTK, ETK, DK, DR, RP, KP motifs
crit1   <- which(sub1a$before2 %in% c("T"))
crit2   <- which(sub1a$before3 %in% c("D", "E"))
crit3   <- which(sub1a$last_second %in% c("T"))
crit4   <- which(sub1a$last_third %in% c("D", "E"))
crit5   <- which(sub1a$first %in% c("P"))
crit6   <- which(sub1a$after1 %in% c("P"))
crit7   <- which(sub1a$before2 %in% c("D"))
crit8   <- which(sub1a$last_second %in% c("D"))
filter1 <- crit1[crit1 %in% crit2] # DTK, ETK
filter2 <- crit3[crit3 %in% crit4] # DEK, ETK

sub2    <- sub1a[which(!1:length(sub1a$peptide) %in% unique(c(filter1, filter2, crit5, crit6, crit7, crit8))),] #
plot3   <- table(sub2$len)
sub3    <- sub2[which(!sub2$before3 %in% c("R", "K") &
                                   !sub2$before2 %in% c("R", "K") &
                                   !sub2$after1  %in% c("R", "K") &
                                   !sub2$after2  %in% c("R", "K")),]
plot4   <- table(sub3$len)

# Visualize bar plot Figure 1B
barplot(plot1, ylim = c(0,65000), col = "red")
barplot(plot2, add = T, col = "orange2")
barplot(plot2a, add = T, col = "orange2")
barplot(plot3, add = T, col = "yellow")
barplot(plot4, add = T, col = "green")



length(unique(sub1a[!is.na(sub1a$V1),"entry"]))

write.csv(combined, "human-trypsinome.csv", row.names = F)
###
files       <- list.files("data/prest_files")
prest_files <- lapply(files, function(x){read.csv(paste0("data/prest_files/", x), header = F)})
prests      <- data.frame(rbindlist(prest_files))

# Proteotypic peptides (fully tryptic, met = 0, unique)

layout(c(1:2))
barplot(c(sum(table(table(sub1a[!is.na(sub1a$V1),"entry"]))[2:45], na.rm = T),
          sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "entry"]))[2:26], na.rm = T)), ylim = c(0,20000), col = "#AFD37E")
barplot(c(sum(table(table(sub1a[!is.na(sub1a$V1),"entry"]))[3:45], na.rm = T),
          sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "entry"]))[3:26])), ylim = c(0,20000), add = T, col = "#98C857")
barplot(c(sum(table(table(sub1a[!is.na(sub1a$V1),"entry"]))[4:45], na.rm = T),
          sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "entry"]))[4:26])), ylim = c(0,20000), add = T, col = "#8CC816")

barplot(c(sum(table(table(sub1a$V1))[2:40], na.rm = T),
          sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "V1"]))[2:40], na.rm = T)))
barplot(c(sum(table(table(sub1a$V1))[3:40], na.rm = T),
          sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "V1"]))[3:40], na.rm = T)), add = T)
barplot(c(sum(table(table(sub1a$V1))[4:40], na.rm = T),
         sum(table(table(sub1a[which(sub1a$V1 %in% prests$V1), "V1"]))[4:40], na.rm = T)), add = T)

