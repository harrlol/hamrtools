library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#arguments: longdf.csv, geneannotation.csv, out dir

args=commandArgs(trailingOnly=TRUE)

g <- fread(args[2], stringsAsFactors = TRUE)
longdf <- fread(args[1], stringsAsFactors = TRUE)

modDistributionSMP <- function(df) {
  t <- df%>%
    filter(lap_type=="gene")
  c <- intersect(unique(t$gene), g$V4)
  toJ <- g%>%filter(V4 %in% c)%>%select(gene=V4, start=V2,end=V3)
  temp <- left_join(t, toJ, by="gene")
  temp%>%
    mutate(norm_pos = (pos-start)/(end-start)*1000)%>%
    ggplot(aes(norm_pos, color=mod))+
    geom_density()+
    facet_wrap(~seq_tech+genotype)+
    theme(plot.caption = element_text(hjust = 0))+
    labs(title="Distribution Landscape of Each Modification Type by Sample Group",
         caption="Mod position is normalized out of 1000, 0 is closer to 5' end of a gene. Empty plot indicates <2 instances present.")+
    xlab("Relative Position")+
    ylab("Frequency")+
    scale_color_manual(values=cbPalette)
}

modDistributionTOT <- function(df) {
  t <- df%>%
    filter(lap_type=="gene")
  c <- intersect(unique(t$gene), g$V4)
  toJ <- g%>%filter(V4 %in% c)%>%select(gene=V4, start=V2,end=V3)
  temp <- left_join(t, toJ, by="gene")
  temp%>%
    mutate(norm_pos = (pos-start)/(end-start)*1000)%>%
    ggplot(aes(norm_pos))+
    geom_density()+
    facet_wrap(~mod)+
    theme(plot.caption = element_text(hjust = 0))+
    labs(title="Distribution Landscape of Each Modification Type Overall",
         caption="Mod position is normalized out of 1000, 0 is closer to 5' end of a gene. Empty plot indicates <2 instances present.")+
    xlab("Relative Position")+
    ylab("Frequency")+
    scale_color_manual(values=cbPalette)
}

modDistributionTOT(longdf)
ggsave(paste0(args[3],"/mod_distribution_along_gene_overall.pdf"), width = 10, height = 8, units = "in")
modDistributionSMP(longdf)
ggsave(paste0(args[3],"/mod_distribution_along_gene_samplegroup.pdf"), width = 10, height = 8, units = "in")
