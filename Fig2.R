library(gwasglue)
#library(devtools)
#devtools :: install_github("mrcieu/gwasglue")
library(dplyr)
library(coloc)
##install.packages("devtools")
devtools::install_github("mglev1n/locusplotr")
library(locusplotr)
library(data.table)
setwd("/home/single1/allbrain/MR")
ctxsa<-fread("ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.txt.gz")
ctxth<-fread("ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.txt.gz")
library(readxl)
eur<-read_excel("eur.xls")
eas<-read.csv("eas.csv")

# Process data to obtain the necessary information for colon analysis
ctxsa$varbeta <- ctxsa$SE^2
ctxth$varbeta <- ctxth$SE^2
eas$sdY = sqrt((eas$beta^2) / (eas$N - 1))
eas$varbeta <- eas$se^2
eur$N<-eur$ncase.exposure
eur$sdY = sqrt((eur$beta.exposure^2) / (eur$N - 1))
eur$varbeta <- eur$se.exposure^2
eur$pval.exposure <- as.numeric(as.character(eur$pval.exposure))
eas$pval <- as.numeric(as.character(eas$pval))

# Create a list and convert it to the format required by the coloc package
ctx1_list <- list(snp = ctxsa$SNP, 
                  beta = ctxsa$BETA1, 
                  varbeta = ctxsa$varbeta, 
                  N = ctxsa$N,
                  MAF=ctxsa$FREQ1,
                  type = "quant") 

ctx2_list <- list(snp = ctxth$SNP, 
                  beta = ctxth$BETA1, 
                  varbeta = ctxth$varbeta, 
                  N = ctxth$N,
                  MAF=ctxth$FREQ1,
                  type = "quant")  

eur_list <- list(snp = eur$SNP, 
                 beta = eur$beta.exposure, 
                 varbeta = eur$varbeta, 
                 N = eur$N,
                 sdY=eur$sdY,
                 type = "quant") 

eas_list <- list(snp = eas$SNP, 
                 beta = eas$beta, 
                 varbeta = eas$varbeta, 
                 N = eas$N,
                 sdY=eas$sdY,
                 type = "quant")  

# Perform colocalization analysis
# Extract the results of colocalization analysis into a new dataframe
# Sort the results by column PP.H4
result1 <- coloc.abf(eur_list,ctx1_list)
result2 <- coloc.abf(eur_list,ctx2_list)
result3 <- coloc.abf(eas_list,ctx1_list)
result4 <- coloc.abf(eas_list,ctx2_list)
res1 <- data.frame(result1$results)
res2 <- data.frame(result2$results)
res3 <- data.frame(result3$results)
res4 <- data.frame(result4$results)
res1 <- res1[order(res1$SNP.PP.H4), ]
res2 <- res2[order(res2$SNP.PP.H4), ]
res3 <- res3[order(res3$SNP.PP.H4), ]
res4 <- res4[order(res4$SNP.PP.H4), ]
# Select the SNP with the largest p-value in PP.H4 and proceed to the next step of plotting
p1<-gg_locusplot(
  df = eur,
  lead_snp = "rs79600142",
  rsid = SNP,
  chrom = chr.exposure,
  pos = pos.exposure,
  ref = effect_allele.exposure,
  alt = other_allele.exposure,
  p_value = pval.exposure,
  plot_distance = 5e+05,
  plot_genes = TRUE,
  population = "EUR",  
  genome_build = "GRCh37"
)
print(p1)

p2<-gg_locusplot(
  df = eur,
  lead_snp = "rs62621252",
  rsid = SNP,
  chrom = chr.exposure,
  pos = pos.exposure,
  ref = effect_allele.exposure,
  alt = other_allele.exposure,
  p_value = pval.exposure,
  plot_distance = 5e+05,
  plot_genes = TRUE,
  population = "EUR", 
  genome_build = "GRCh37"
)
print(p2)

p3<-gg_locusplot(
  df = eas,
  lead_snp = "rs16940665",
  rsid = SNP,
  chrom = chr,
  pos = pos,
  ref = major,
  alt = other_allele,
  p_value = pval,
  plot_distance = 5e+05,
  plot_genes = TRUE,
  population = "EAS", 
  genome_build = "GRCh37"
)
print(p3)

p4<-gg_locusplot(
  df = eas,
  lead_snp = "rs12373139",
  rsid = SNP,
  chrom = chr,
  pos = pos,
  ref = major,
  alt = other_allele,
  p_value = pval,
  plot_distance = 5e+05,
  plot_genes = TRUE,
  population = "EAS",  
  genome_build = "GRCh37"
)
print(p4)
