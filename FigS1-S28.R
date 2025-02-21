library(TwoSampleMR)
library(MRPRESSO)
library(data.table)
setwd("/path/")

# Transform GWAS data of PD in European populations 
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")
# install.packages("devtools")
# devtools::install_github( "MRCIEU/gwasglue" )
library(VariantAnnotation)
library(gwasglue)
EUR_exposure = readVcf("ieu-b-7.vcf.gz")
EUR_exposure = gwasvcf_to_TwoSampleMR(vcf = EUR_exposure,type = "exposure")
head(EUR_exposure)
write.csv(EUR_exposure, file="PD_EUR_exposure.csv",sep="\t",quote = F,row.names = F)

# Export the exposure file to the TwoSampleMR package under R file 
# Clumping process
setwd("path/R/win-library/4.3/TwoSampleMR")
exp_dat<-system.file("PD_EUR_exposure.csv",package = "TwoSampleMR")

# Clumping setting
exp_dat_clumped<-read_exposure_data(
  exp_dat,sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  clump = FALSE
)
exp_dat_clumped<-clump_data(
  exp_dat_clumped,
  clump_kb=10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 2,
  pop = "EUR"
)

# Filter exp SNPs 
exp_dat_clumped<-subset(exp_dat_clumped,pval.exposure<5e-08)
# read outcome data
out_dat<-fread("ENIGMA3_mixed_se_wo_Mean_bankssts_surfavg_20190429.txt.gz")
merge<-merge(exp_dat_clumped,out_dat,by.x="SNP",by.y="SNP")
write.csv(merge,"merge_se_wo_Mean_xxx_SurfArea.csv") # rename different region files
out_dat<-read_outcome_data(
  snps = merge$SNP,
  filename = "merge_se_wo_Mean_xxx_SurfArea.csv", # rename different region files
  sep=",",
  snp_col = "SNP",
  beta_col = "BETA1",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FREQ1",
  pval_col = "P"
)
# Filter out SNPs 
out_dat<-subset(out_dat,pval.outcome>0.05)
# Use the harmonise function to remove the palindromic sequence
dat <- harmonise_data(exposure_dat = exp_dat_clumped, outcome_dat = out_dat, action = 2)

# Perform MR
result <-generate_odds_ratios(mr(dat))
write.table(result, file="MR-result.xls",sep="\t",quote=F)
res_single <- mr_singlesnp(dat)
singlesnpOR <-generate_odds_ratios(res_single)
write.table(singlesnpOR, file="singlesnpOR.xls",sep="\t",quote=F)

# scatter plot
mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode")),dat)
# funnel plot
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
# forest plot
mr_forest_plot(singlesnp_results = mr_singlesnp(dat))
# leave one out plot
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))

# Pleiotropic test
pleiotropy <- mr_pleiotropy_test(dat)
write.table(pleiotropy, file="pleiotropy-test.xls",sep="\t",quote=F)
# heterogeneity test
heterogeneity<-mr_heterogeneity(dat)
write.table(heterogeneity, file="Q-test.xls",sep="\t",quote=F)
# mr_presso detects outliers
mr_presso<-run_mr_presso(dat,NbDistribution = 10000)
write.table(mr_presso, file="outliers.xls",sep="\t",quote=F)
