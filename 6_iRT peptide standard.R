library(stringr)
library(ggplot2)
library(reshape)

df      <- read.csv("data/20160705_rt_peptides_2.csv")
filter  <- c("GSHMASLAEAK", "SQTPAEDTVK", "YGVSDYHKNLINNAK", "ELDKYGVSDYHKNLINNAK")
df      <- df[which(!df$Peptide %in% filter),] # delete peptides*
df$grad <- unlist(lapply(df$Replicate, function(x){str_split(x, "_")[[1]]}[2]))
aggregate(df$Peptide.Retention.Time, list(df$Peptide, df$grad), median)
rt_cal <- aggregate(df$Peptide.Retention.Time, list(df$Peptide, df$grad), median)
rt_cal <- rt_cal[rt_cal$Group.2 == "90min",]

irt <- function(retention_times){
  return((retention_times - min(retention_times))/(max(retention_times) - min(retention_times))*100)
}

rt_cal$irt <- irt(rt_cal$x)
plot <- merge(df, rt_cal[,c("Group.1", "irt")], by.x = c("Peptide"), by.y = c("Group.1"), all.x = T, all.y = F)

cv <- function(x){
  return(round(sd(x)/as.numeric(mean(x))*100,2))
}

ggplot(plot, aes(y = as.numeric(Peptide.Retention.Time), x = as.numeric(irt), color = factor(Replicate))) + geom_point(position = ) + theme +
  xlab("iRT") + ylab("Retention Time (min)") + geom_smooth(method = lm)

grad30 <- plot[which(plot$grad == "30min"),]
grad30 <- grad30[!duplicated(grad30[,c("Peptide", "Replicate")]),]
temp <- data.frame(name = unique(grad30$Peptide), a = rstandard(lm(grad30[which(grad30$Replicate == "FE_30min_2ul_20160705_r1"),c("Peptide.Retention.Time", "irt")])),
           b = rstandard(lm(grad30[which(grad30$Replicate == "FE_30min_2ul_20160705_r2"),c("Peptide.Retention.Time", "irt")])),
           c = rstandard(lm(grad30[which(grad30$Replicate == "FE_30min_2ul_20160705_r3"),c("Peptide.Retention.Time", "irt")])))

a       <- merge(temp, standard, by.x = "name", by.y = "Peptide", all.x = T, all.y = F)
a$a_irt <- a[,2] + a[,5]
a$b_irt <- a[,3] + a[,5]
a$c_irt <- a[,4] + a[,5]
a$cv    <- apply(a[,6:8],1, cv)

write.csv(a, "output/30min.csv")

rep <-  unique(plot$Replicate)
r2 <- list()
for( i in 1: length(rep))(
  r2[[i]] <- summary(lm(plot[which(plot$Replicate %in% rep[i]), c("Peptide.Retention.Time", "irt")]))$r.squared
)

theme <- theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               panel.background = element_blank(), 
               legend.title=element_blank(),
               axis.line = element_line(color = "black"), 
               axis.ticks.x = element_blank(),
               axis.title = element_text(),
               plot.title = element_text(size = 24),
               axis.ticks.y = element_blank())