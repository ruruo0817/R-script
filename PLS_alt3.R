setwd("E:/Rpractice/meta/tem_article/PLS/alt3_pls/diversity")
# Install the plspm package
# install.packages('devtools')
# devtools::install_github('gastonstat/plspm')
# For detailed results, refer to the user manuals of the plspm package:
# Full manual, 235 pages: https://www.gastonsanchez.com/PLS_Path_Modeling_with_R.pdf 
# Brief manual, 10 pages: https://rdrr.io/cran/plspm/f/inst/doc/plspm_introduction.pdf 
# View the path coefficient parameter estimates and related statistical information
# Load the plspm package
library(plspm)
library(ggplot2)
# I. Aquatic structure equation model
# Read data
# Read data
phy <- read.delim('water903_slope.txt', sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
head(phy)
dat <- subset(phy, type == "water") # Specify latent variables, store the relationship between variables and latent variables in R using a list (list)
# First model
dat_blocks <- list(
  climate = c("PRE", "TEM"),
  NPP = "NPP", 
  # NDVI = "NDVI",
  local = c('WT', 'DO', 'pH', "TN", "TP", "NH4_N"), 
  community = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
                "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes",
                "Others", "Planctomycetota",
              "Proteobacteria", "Verrucomicrobiota")
# diversity = "shannon_diversity"
)
dat_blocks
# Describe the relationships between latent variables through a 0-1 matrix, where 0 means no relationship between variables, and 1 means there is a relationship

climate <- c(0, 0, 0, 0)
NPP <- c(0, 0, 0, 0)
# NDVI <- c(0, 0, 0, 0, 0)
local <- c(1, 1, 0, 0)
community <- c(1, 1, 1, 0)
# community <- c(1, 1, 1, 1, 0)

dat_path <- rbind(climate, NPP, local, community)
colnames(dat_path) <- rownames(dat_path)
dat_path
colnames(dat_path) <- rownames(dat_path)
dat_modes <- rep('A', 4)

# Second model
dat_blocks <- list(
 NPP = "NPP", 
 NDVI = "NDVI",
 climate = c("PRE", "TEM"),
 local = c('WT', 'DO', 'pH', "TN", "TP", "NH4_N"), 
 community = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota",
               "Chloroflexi", "Cyanobacteria", "Firmicutes", "Gemmatimonadetes" 
               , "Nitrospirae", "Others", "Patescibacteria", "Planctomycetota",
               "Proteobacteria", "Verrucomicrobiota"),
diversity = "shannon_diversity"
)
dat_blocks
# Describe the relationships between latent variables through a 0-1 matrix, where 0 means no relationship between variables, and 1 means there is a relationship
NPP <- c(0, 0, 0, 0, 0, 0)
NDVI <- c(0, 0, 0, 0, 0, 0)
climate <- c(0, 0, 0, 0, 0, 0)
local <- c(1, 1, 1, 0, 0, 0)
community <- c(1, 1, 1, 1, 0, 0)
diversity <- c(1, 1, 1, 1, 1, 0)
dat_path <- rbind(NPP, NDVI, climate, local, community, diversity)
colnames(dat_path) <- rownames(dat_path)
dat_path
# Specify the causal relationship,可选择 A (representing columns are the cause of rows) or B (representing rows are the cause of columns)
dat_modes <- rep('A', 6)
dat_modes

## A simple PLS-PM, more parameter details ?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)
capture.output(summary(dat_pls), file = "water_pls_summary.txt")
dat_pls$path_coefs
dat_pls$inner_model

# Output path coefficients and statistical information
# Create an empty data frame
result <- data.frame()
# Add the results of $community and $diversity to the data frame one by one

result <- rbind(result, dat_pls$inner_model$community)
write.csv(result, "water_results.csv", row.names = TRUE)
# The goodness-of-fit value can help assess model fit
gof <- round(dat_pls$gof, 3)
# View the saved path diagram as a JPEG file
jpeg("inpath_water.jpeg", width = 800, height = 600, res = 120)
# Plot the path diagram
innerplot(dat_pls, colpos = "orange", colneg = 'skyblue', show.values = TRUE, lcol = 'black', box.lwd = 0)
# Add title
title(main = paste("water PLS-PM (gof =", gof, ")"))
dev.off() # Close the graphics device
# View the status of external and internal latent variables
dat_pls$inner_summary
# View the impact status of variables
dat_pls$effects
# View the relationship between observed variables and latent variables
dat_pls$outer_model
# View the saved path diagram as a JPEG file
jpeg("outpath_water.jpeg", width = 1000, height = 700, res = 150)
# You can use outerplot() to plot a structure similar to a path diagram, details ?outerplot
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'orange', colneg = 'skyblue', show.values = TRUE, lcol = 'black')
# Add title weights
title(main = paste("water outpath"))
dev.off() # Close the graphics device

##____2.for loop
# Create a vector containing all group names
groups <- c( "L", "M", "H")
# Use a for loop to iterate through all groups, and perform structural equation modeling for each category
# Define group names
unique(dat$alt3)
groups <- unique(dat$alt3)
for (i in groups) {
  # Filter data based on group name
  dat_i <- dat[grepl(i, dat$alt3), ]
  # Run structural equation model
  pls_i <- plspm(dat_i, dat_path, dat_blocks, modes = dat_modes)
  summary(pls_i)
  capture.output(summary(pls_i), file = paste0("water_pls.sammary_", i, ".txt"))
  pls_i$path_coefs
  pls_i$inner_model
  
  # Create an empty data frame
  result <- data.frame()
  # Add the results of $community and $diversity to the data frame one by one
  result <- rbind(result, pls_i$inner_model$community)
  result <- rbind(result, pls_i$inner_model$diversity)
  filename <- paste0("water_results", i, ".csv")
  write.csv(result, filename, row.names = TRUE)
  
  # Output path diagram
  jpeg(paste("inpath_water_", i, ".jpeg"), width = 800, height = 600, res = 120)
  innerplot(pls_i, colpos = "orange", colneg = 'skyblue', show.values = TRUE, lcol = 'black', box.lwd = 0)
  title(main = paste(i, "-water_PLS(gof =", round(pls_i$gof, 2), ")"))
  dev.off()
  pls_i$inner_summary
  # View the impact status of variables
  pls_i$effects
  # View the relationship between observed variables and latent variables
  pls_i$outer_model
  # View the saved path diagram as a JPEG file
  # Output path diagram
  jpeg(paste( "outpath_water_", i, ".jpeg"), width = 800, height = 