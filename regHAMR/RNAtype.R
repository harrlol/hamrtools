library(ggplot2)
library(dplyr)
library(data.table)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args=commandArgs(trailingOnly=TRUE)

df <- fread(args[1])
dir <- dirname(args[1])

a <- unique(df$genotype)
b <- unique(df$seq_tech)
g <- expand.grid(a,b)

for (i in (1:nrow(g))) {
  # Create Data
  data <- df%>%
    filter(genotype == g[i,1] & seq_tech == g[i,2] & lap_type %in% c("ncRNA", "gene"))%>%
    group_by(bio)%>%
    summarize(count=n())
  
  # Compute the position of labels
  data <- data %>% 
    arrange(desc(bio)) %>%
    mutate(prop = round(count / sum(data$count) *100, digits=1))%>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  ggplot(data, aes(x="", y=prop, fill=bio)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=cbPalette)+
    theme_void()+
    labs(title=paste("Mod Distribution in RNA Subtypes in", g[i,1], g[i,2], sep = " "))
  
  ggsave(paste(dir,"/RNA_type_", g[i,1], "_", g[i,2], ".pdf", sep=""), width = 10, height = 8, units = "in")
}
 