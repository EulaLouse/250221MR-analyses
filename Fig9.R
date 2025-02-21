rm(list=ls())
# Set up work environment
setwd("/path/")


# install.packages("ggplot2")
# install.packages("ggalluvial")
library(ggplot2)
library(ggalluvial)

# read data
data <- read.csv("Sankey_dat.csv",header = T, check.names = F)

# transfer format
df <- to_lodes_form(data[,1:ncol(data)],
                    axes = 1:ncol(data),
                    id = "value")
print(df)

# Sankey diagram
col<- rep(c('#7b4173', '#8ca252', '#cedb9c', '#ad494a',
            '#9e9ac8', '#d6616b', '#e7cb94','#598c14', 
            '#bcbddc', '#393b79', '#637939', '#843c39',
            '#de9ed6','#9ecae1','#e7969c','#e7ba52','#8c6d31','#5254a3',
            '#b5cf6b','#bd9e39','#dadaeb','#3182bd',
            '#5254a3','#a55194','#c6dbef','#ce6dbd',
            '#9c9ede','#6b6ecf','#756bb1','#6baed6'), 3)

pdf("test.pdf",width = 8, height = 6)
ggplot(df, aes(x = x, fill=stratum, label=stratum,
               stratum = stratum, alluvium  = value))+
  geom_flow(width = 0.2,
            curve_type = "sigmoid",
            alpha = 0.5,
            color = 'white',
            size = 0.2)+
  geom_stratum(width = 0.28)+
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+
  theme_void()+
  theme(axis.line=element_line(linetype=1,color="grey",size=1.5),
        plot.title = element_text(size=15,hjust = 0.5), 
        legend.position = 'none')
dev.off()