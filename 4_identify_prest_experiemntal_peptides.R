library(data.table)
library(stringr)
prestpep <- list.files("data/prest_files/")
files <- list.files("data/prest_peptides/")
all_pep <- lapply(files, function(x){readRDS(paste0("data/prest_peptides/",x))})
all_pep <- rbindlist(all_pep)
all_pep$Start.position <- as.numeric(all_pep$Start.position)

all_pep$y_ions <- unlist(lapply(all_pep$Matches, function(x){str_count(x, "y")}))
all_pep$b_ions <- unlist(lapply(all_pep$Matches, function(x){str_count(x, "b")}))
all_pep$fragments <- (all_pep$y_ions + all_pep$b_ions)
hist(all_pep$fragments)


