#install.packages("pheatmap")
library(pheatmap)
inputFile="heatmapdata.txt"
groupFile="heatmapgroup.txt"
setwd("/path/")      
#input data file
dat=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
#input sample group file
group=read.table(groupFile,header=T,sep="\t",row.names=1,check.names=F)    

pheatmap(dat,
         annotation=group,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("#b10b26","#fdbb6d" ,"#ebf7e3", "#72abd0","#313695"))(50),
         show_colnames = T,
         scale="none",  
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)

