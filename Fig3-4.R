library(ggplot2)
# install.packages("devtools")
# devtools::install_github("adayim/forestploter") 
library(forestploter)
# Set up work environment
setwd("path/Fig2-3.R")

pacman::p_load(forestploter,readxl,haven)
forest<- read_excel("forest_SA_ctx_region.xlsx")
forest$Subgroups <-ifelse(is.na(forest$"Event"),
                          forest$Subgroups,
                          paste0(" ", forest$Subgroups))
#Create an empty column to store the graphic part of the forest map, with a length of 20 
forest$` ` <- paste(rep(" ", 20), collapse ="")
# NA needs to be converted to blank
forest$'OR(95%CI)AUTO' <- ifelse(is.na(forest$OR), "",sprintf("%.2f (%.2f to %.2f)",forest$OR, forest$LowerCI, forest$UpperCI))
forest$"Event" <- ifelse(is.na(forest$"Event"), "", forest$"Event")
forest$"OR(95%CI)" <- ifelse(is.na(forest$"OR(95%CI)"), "", forest$"OR(95%CI)")
plot1 <- forest(forest[, c(1:3,7,8)],
                est = forest$OR,
                lower=forest$LowerCI,
                upper=forest$UpperCI,
                ci_column=5,
                ref_line=1,
                xlim = c(-20,15),ticks_at = c(-20,-15,-10,-5,0,5,10,15))
plot1
tm <- forest_theme(base_size = 12, base_family = "",
                   ci_pch = 16, ci_col = "#4575b4", ci_lty = 1, ci_lwd =1.5,ci_Theight=0.2,
                   legend_name = "Group", legend_position = "right", legend_value="",xaxis_lwd= 0.6, xaxis_cex= 1,
                   footnote_cex = 0.6, footnote_fontface = "plain", footnote_col="black")
# change theme
plot2 <- forest(forest[, c(1:3,7,8)],
                est = forest$OR,
                lower = forest$LowerCI, upper = forest$UpperCI, ci_column =5,
                ref_line = 1,
                xlim = c(-20,15),ticks_at = c(-20,-15,-10,-5,0,5,10,15), 
                theme = tm)
plot2


# Draw multi group forest plot
pacman::p_load(forestploter,readxl,haven)
forest<- read_excel("forest_SA_multi_region.xlsx")
# Empty column "TotalY1" to store other information in the panel
forest$Regions <- ifelse(is.na(forest$"TotalY1"),
                         forest$Regions,
                         paste0(" ", forest$Regions))
# NA needs to be converted to blank
forest$TotalY1<-ifelse(is.na(forest$TotalY1), "", forest$TotalY1)
#Create an empty column to store the graphic part of the forest plot, with a length of 30
forest$`n1`<- ifelse(is.na(forest$TotalY1), "", forest$TotalY1)
forest$' ' <- paste(rep(" ", 30), collapse="")
#Create another two empty columns to store the results of different panels, with a length of 20
forest$"Beta1(95%CI)" <- ifelse(is.na(forest$"Beta1(95%CI)"), "", forest$"Beta1(95%CI)")
forest$"Beta2(95%CI)" <- ifelse(is.na(forest$"Beta2(95%CI)"), "", forest$"Beta2(95%CI)")
forest$interaction_p <- c(0.858,0.594,0.763)
tm <- forest_theme(base_size =10, base_family ="",
                   ci_pch =c(16,18), ci_col =c("#e7969c","#81b7d7"),
                   ci_lty = 1,ci_lwd = 1.5, ci_Theight =0.4,
                   legend_name = "Group", legend_position = "right",
                   legend_value = c("European", "East Asian"), 
                   xaxis_lwd = 0.6, xaxis_cex=1,
                   refline_lwd = 1, refline_lty = "dashed", refline_col="#4574b4",
                   vertline_lwd = 1,
                   vertline_lty = c("dashed","dotted"),
                   vertline_col = c("#d6604d","#bababa"),
                   summmary_fill="#4574b4",summary_col = "#4575b4",
                   footnote_cex=0.6,footnote_fontface="plain",footnote_col="black")
plot4<-forest(forest[,c(1:3,14,11,12,15)],
              est=list(forest$OR,forest$OR2),
              lower=list(forest$LowerCI,forest$LowerCI2),
              upper=list(forest$UpperCI,forest$UpperCI2),
              ci_column=c(4),
              ref_line=1,
              xlim=c(0,10),
              ticks_at=c(0,2,4,6,8,10),
              add_interaction_p = TRUE,  # This parameter is used to add the interaction p-value column
              interaction_p_col = "interaction_p",
              theme=tm)
plot4


# calculate interaction p
library(dplyr)
data <- data.frame(
  group = factor(c("Group 1", "Group 2")),  # group
  OR = c(1.0054, 0.995),  # Odds Ratios
  SE = c(0.0063, 0.0019)  # Standard Errors
)

# Compute the logarithmic OR
data <- data %>%
  mutate(
    logOR = log(OR),
    var = SE^2  # Calculate variance
  )

# Calculating standard error for interaction effects
# The standard error of interaction is equal to the square root of the sum of the variances of the ancestry-specific effect estimates
interaction_se <- sqrt(sum(data$var))

# Calculate interaction effects (on logarithmic scale)
interaction_effect <- data$logOR[1] - data$logOR[2]  # "logOR of Group1" minus "logOR of Group2"

# Calculate the z-value of the interaction
z_value <- interaction_effect / interaction_se

# p-values for interaction were calculated using two-tailed tests
interaction_p_value <- 2 * (1 - pnorm(abs(z_value)))  

# print p
print(paste("interaction-p:", interaction_p_value))

