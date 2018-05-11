# Denna som används för publikationen

library(data.table)
library(stringr)
combined <- readRDS("data/20180503_peptides_all_prest_combined_flanking.RDa")
prest <- list.files("data/prest_files/")
prests <- lapply(prest, function(x){read.csv(paste0("data/prest_files/",x), header = F)})
all_prests <- rbindlist(prests)
files <- list.files("data/prest_peptides/")
all_pep <- lapply(files, function(x){readRDS(paste0("data/prest_peptides/",x))})
all_pep <- rbindlist(all_pep)
all_pep$Start.position <- as.numeric(all_pep$Start.position)

all_pep$y_ions <- unlist(lapply(all_pep$Matches, function(x){str_count(x, "y")}))
all_pep$b_ions <- unlist(lapply(all_pep$Matches, function(x){str_count(x, "b")}))
all_pep$fragments <- (all_pep$y_ions + all_pep$b_ions)
all_pep$hyd <- unlist(lapply(all_pep$Sequence, function(x){hydrophobicity(x, scale = "KyteDoolittle")}))
all_pep$Sequence <- paste(all_pep$Sequence)
combined$experiemental <- combined$peptide %in% all_pep$Sequence

combined$hyd <- unlist(lapply(combined$peptide, function(x){hydrophobicity(x, scale = "KyteDoolittle")}))
hist(combined$hyd[which(combined$experimental == T)])



quanto <- combined[which(combined$prot_mod == "" & combined$prot_pep == "" & combined$met == 0),] # subset to get quantotypic
proto <- combined[combined$met == 0,]
a <- all_pep[all_pep$Sequence %in% quanto$peptide,]
b <- all_pep[all_pep$Sequence %in% proto$peptide,]
e <- all_pep[all_pep$Sequence %in% combined$peptide,]
length(unique(a$Proteins))
length(unique(b$Proteins))
length(unique(e$Proteins))
# 26840 prests in total
# 12849 yield quantotypic peptides
# 14186 yield proteotypic peptides
# 15682 fully tryptic
# 21583 yield any unique peptide 

# map to chromosomes

write.csv(b, "proteotypic-peptides.csv", row.names = F)

df <- read.table("data/hprr-by-chromosome.xls", header = T)
prot <- data.frame(table(b[!duplicated(b[,c("Sequence", "Leading.razor.protein")]), c("Leading.razor.protein")]))

df2 <- merge(df, prot, by.x = "PrEST", by.y = "Var1", all = T)

res <- aggregate(df2$Freq, list(df2$Gene, df2$Chromosome), sum, na.rm = T)
res$temp <- unlist(lapply(res$Group.2, function(x){unlist(str_split(x, ";"))[1]}))
res$temp <- factor(res$temp, levels = c('1','2','3','4','5','6','7','8','9',
                                        '10','11','12','13','14','15','16','17','18','19',
                                        '20','21','22', 'X','Y', "MT", 'Unmapped'))
barplot(table(res[, "temp"]) / table(res[, "temp"]) * 100, col = "#FFFFFF", ylim = c(0,100))                                           
barplot(table(res[which(res$x > 0), "temp"]) / table(res[, "temp"]) * 100, col = "#AFD37E", add = T)
barplot(table(res[which(res$x > 1), "temp"]) / table(res[, "temp"]) * 100, add = T, col ="#98C857")
barplot(table(res[which(res$x > 2), "temp"]) / table(res[, "temp"]) * 100, add = T, col = "#8CC816")

table(df$Chromosome)

table(res$temp)


write.csv(, "proteotypic.csv")

pie(c(26840-21583,21583-15682,15682-14186,14186))

26840 - 21583
21583 - 15682
15682 - 14186
14186 - 12849
sum(table(table(combined[which(combined$experiemental == T & combined$met == 0),]$entry))[4:20], na.rm = T)
sum(table(table(combined[which(combined$experiemental == T & combined$met == 0),]$V1))[2:20], na.rm = T)
head(combined)
sum(table(table(combined[which(!is.na(combined$V1) & combined$met == 0),]$V1))[2:20], na.rm = T)
mc1 <- unique(all_pep[all_pep$Missed.cleavages == 1]$Sequence)

write.csv(mc1, "missed-cleaved-peptides.csv")

length(grep("RR|KK|RK|KR", mc1))
length(grep("DTR|DTK|DK|DR|ETR|ETK|ER|EK", mc1))
length(grep("K.K$|R.K$|R.R$|K.R$", mc1))
length(grep("RP|KP", mc1))
length(mc1)
# length(mc1)
pie(c(1224, 3502, 8999, 3191, (28365-3502-8999-3191-1224)))
hist(all_pep$iRT)          


# plot figure 2B
uni <- unique(b$Sequence)
barplot(table(nchar(uni)), xlab = "Number of amino acids", ylab = "Number of PrESTs")

# plot figure 2C

b <- data.frame(b)
uni2 <- unique(b[, c("Sequence", "Proteins")])
bplot <- table(table(uni2$Proteins))
barplot(bplot, ylim = c(0,8000))
text(x = 1:9, y = bplot, label = bplot, pos = 3, cex = 0.8, col = "red")
sum(bplot)

proto <- data.frame(proto)
bplot <- table(table((proto[which(proto$experiemental == T), "entry"])))[2:16]
barplot(bplot, ylim = c(0,4000))
text(x = 1:9, y = bplot, label = bplot, pos = 3, cex = 0.8, col = "red")
sum(bplot)

hist(b$fragments, xlim =c(0,70))
