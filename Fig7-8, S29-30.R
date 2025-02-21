#BiocManager::install("clusterProfiler")
#BiocManager::install("topGO")
#BiocManager::install("Rgraphviz")
#BiocManager::install("pathview")
#BiocManager::install("org.Hs.eg.db")
#install.packages("openxlsx")
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(openxlsx)
library(readxl)
# Set up work environment
setwd("/path/")
data<-read_xlsx("gene.xlsx",sheet=1)
head(data)
keytypes(org.Hs.eg.db)
#Convert gene ID
entrezid_all<-mapIds(
  x=org.Hs.eg.db,
  keys=data$gene_name,
  keytype="SYMBOL",
  column="ENTREZID"
)
#Remove the NA value in the conversion result
entrezid_all<-na.omit(entrezid_all)
#Convert results to data frame format
entrezid_all<-data.frame(entrezid_all)
head(entrezid_all)

# GO enrichment analysis
GO_enrich=enrichGO(gene = entrezid_all[,1],
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = T)
GO_enrich=data.frame(GO_enrich)
# export results
write.csv(GO_enrich,file="GOpathways.csv")

# KEGG enrichment analysis
KEGG_enrich=enrichKEGG(gene=entrezid_all[,1],
                       keyType = "kegg",
                       pAdjustMethod = 'fdr',
                       organism = "human",
                       qvalueCutoff = 1,
                       pvalueCutoff = 1)
KEGG_enrich=data.frame(KEGG_enrich)
# export results
write.csv(KEGG_enrich,file = "KEGGpathways.csv")

# Wiki enrichment analysis
WIKI_enrich <- enrichWP(gene=entrezid_all[,1], 
                      organism = "Homo sapiens",
                      qvalueCutoff = 1,
                      pvalueCutoff = 1)
write.csv(WIKI_enrich,file = "Wikipathways.csv")

# read the filtered GO pathways file
# GO pathways bar plot
go_enrich<-read.csv("GOpathways.csv",stringsAsFactors = FALSE,fileEncoding = 'GBK')
go_enrich$term<-paste(go_enrich$ID,go_enrich$Description,sep = ":")
go_enrich$term<-factor(go_enrich$term,levels = go_enrich$term,ordered = TRUE)
library(dplyr)
go_enrich_filtered<-go_enrich %>%
  arrange(pvalue)%>%
  group_by(ONTOLOGY)%>%
  slice_head(n=10)%>%
  ungroup()
go_enrich_filtered$term<-factor(go_enrich_filtered$term,levels = rev(go_enrich_filtered$term))
ggplot(go_enrich,aes(x=term,y=EnrichmentScore,fill=ONTOLOGY))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = c("#eb6837","#269d6b","#373a79"))+
  facet_wrap(~ONTOLOGY,ncol = 1,
             scales = "free_y")+
  coord_flip()+
  xlab("GO term")+
  ylab("Enrichment Score")+
  labs(title="GO Results of Three Ontologies")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8,angle = 90,hjust = 1))
go_enrich<-go_enrich[order(go_enrich$pvalue),]

# read the filtered KEGG pathways file
# KEGG pathways bar plot
kk_result<-read.csv("KEGGpathways.csv",stringsAsFactors = FALSE,fileEncoding = 'GBK')
kk_result<-kk_result[order(kk_result$pvalue),]
display_number<-30
kk<-head(kk_result,n=display_number)
kk$order<-factor(rev(as.integer(rownames(kk))),labels=rev(kk$Description))
ggplot(kk,aes(y=order,x=Count,fill=pvalue))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_gradient(low = "red",high = "blue")+
  labs(title = "KEGG Pathways Enrichment",
       x="Gene number",
       y="Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

# KEGG bubble plot
ggplot(kk,aes(y=order,x=EnrichmentScore))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient2(
    low = c("#b10b26","#fdbb6d"),
    mid = "#fffdd1",
    high = c("#72abd0","#313695"),
    midpoint = 0.125,
    name = "pvalue")+  
  scale_size(range = c(3,4))+
  labs(color=expression(Pvalue,size="Count"),
       x="EnrichmentScore(-log10(pvalue))",y="Pathways",title = "KEGG Pathways Enrichment")+
  xlim(0,2.0)+
  theme_bw()

# read the filtered WIKI pathways file
# WIKI pathways bubble plot
wiki_result<-read.csv("Wikipathways.csv",stringsAsFactors = FALSE,fileEncoding = 'GBK')
wiki_result<-wiki_result[order(wiki_result$pvalue),]
display_number<-30
wp<-head(wiki_result,n=display_number)
wp$order<-factor(rev(as.integer(rownames(wp))),labels=rev(wp$Pathway))

ggplot(wp,aes(y=order,x=Enrichment_score))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient2(
    low = c("#b10b26","#fdbb6d"),
    mid = "#fffdd1",
    high = c("#72abd0","#313695"),
    midpoint = 0.025,)+  
  scale_size(range = c(3,4))+
  labs(color=expression(Pvalue,size="Count"),
       x="EnrichmentScore",y="Pathways",title = "Wiki Pathways Enrichment")+
  theme_bw()