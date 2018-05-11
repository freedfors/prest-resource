library(gplots)
library(Metrics)
library(ggplot2)
library(grid)
library(gridExtra)
library(ClassDiscovery)
library(gplots)
library(reshape2)
library(PerformanceAnalytics)
library(Rtsne)
library(scatterplot3d)
library(ggrepel)
library(openxlsx)
library(functional)
library(RColorBrewer)
library(mvtboost)

fill.na <- function(dataframe, columns, cutoff = 0, fill = NA){
  subset <- dataframe[,columns]
  subset[subset <= cutoff] <- fill
  dataframe[,columns] <- subset
  return(dataframe)
}


df <- read.csv(file.choose())
df <- df[-which(df$unifier == "APOE - QPrEST 2" | df$unifier == "C8G - QPrEST 1"),]


df2 <- df[,c("unifier", "sample", "conc.fmol.ul")]
df3 <- dcast(df2, unifier ~ sample , median)

id <- df3[,1]
df3 <- df3[,-1]

# median normalized for correlation analysis
df4 <- df3/apply(df3, 1, median, na.rm  =T)

df3 <- log2(df3)
df3[is.na(df3)] <- 0
df4[is.na(df4)] <- 0

pca         <- prcomp(df4)
rot         <- pca$r
scores      <- pca$x
theme       <- theme_minimal() + theme(panel.grid = element_blank(),
                                       panel.background = element_rect(fill = "white", colour = "grey50"),
                                       axis.text.x = element_text(angle = 90),
                                       legend.title= element_blank(),
                                       legend.position = "none")
sex <- df[!duplicated(df[,1,6]), c(1,6,7)]
sex$id <- gsub("_.*", "", sex$sample.id)
sex <- merge(sex, data.frame(unique(sex$id), rainbow(17)), by.x = "id", by.y = "unique.sex.id.", all = T)
colnames(sex) <- c("id", "sample", "sex", "sample.id", "col")
sex$visit <- gsub(".*_", "", sex$sample.id)
sex <- sex[order(sex$sample),]
sex$sex.col <- "#E66631"
sex[sex$Sex == "m",]$sex.col <- "#A5ADB8"

p1 <- ggplot(data.frame(rot), aes(x = PC1, y = PC2, col = factor(sex[,3]))) + 
  geom_point() + theme 

# plot names instead
#p1 <- ggplot(data.frame(rot), aes(x = PC1, y = PC2, label = colnames(df3), col = factor(sex[,2]))) + 
#  geom_text() + theme 

p1b <- ggplot(data.frame(scores), aes(x = PC1, y = PC2, label = id)) + geom_text() + theme

p2 <- ggplot(data.frame(rot), aes(x = PC1, y = PC3, col = factor(sex[,3]))) + 
  geom_point() + theme

p2b <- ggplot(data.frame(scores), aes(x = PC1, y = PC3, label = id)) + geom_text() + theme

p3 <- ggplot(data.frame(rot), aes(x = PC2, y = PC3, col = factor(sex[,4]))) +
  geom_point() + theme 
p3b <- ggplot(data.frame(scores), aes(x = PC2, y = PC3, label = id)) + geom_text() + theme


# TSNE
tsne_out = Rtsne(t(df3), perplexity = 10)
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])

# print file

pdf("results.pdf", height = 9, width = 3)
grid.arrange(p1, p2, p3)
grid.arrange(p1b, p2b, p3b)
dev.off()

ggplot(tsne_plot, aes(x = x, y = y, col = factor(sex[,1]))) + geom_point() + theme + xlab("tsne.x") + ylab("tsne.y")

# histogram
df$log2conc <- log2(df$conc.fmol.ul)

ord <- aggregate(df$conc.fmol.ul, list(df$unifier), median, na.rm = T)
df$unifier <- factor(df$unifier, levels = ord[order(ord[,2], decreasing = T),1])

# BOXPLOT
ggplot(data.frame(df), aes(x = unifier, y = log10(conc.fmol.ul), col = Sex)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_dodge(width=0.75)) +
  theme +
  xlab("Gene ID")

heatmap.2(cor(df4, method = "spearman", use = "pairwise.complete.obs"),  labCol = sex[,4], labRow = sex[,3], 
          trace = "none", colCol = paste0(sex$rainbow.17.))


heatmap.2(t(as.matrix(df3)),  labRow = sex[,4], labCol = levels(df2$unifier), 
          trace = "none", colRow = paste0(sex$sex.col))
