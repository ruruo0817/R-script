Here is the English translation of the provided R code:

```r
setwd("E:/Rpractice/meta/tem_article/rda")
#_______________Standardize environmental factor data and perform normality test————————
# Read data
phy <- read.delim('phy_env1252.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

head(env) # Read data
envnum <- env[,-3]
head(envnum)
# Transform data, for example, take the logarithm
trans_env <- log(envnum+1)
# Save the standardized results
write.csv(trans_env, "env_selected_transform.csv")
# Perform Shapiro-Wilk normality test
shapiro_results <- lapply(trans_env, shapiro.test)
print(shapiro_results)
# Save the test results as a data frame
shapiro_df <- data.frame(
  Variable = names(shapiro_results),
  p_value = sapply(shapiro_results, function(x) x$p.value),
  W_statistic = sapply(shapiro_results, function(x) x$statistic)
)
# Save the results as a CSV file
write.csv(shapiro_df, "Shapiro-Wilk Test Results.csv", row.names = FALSE)

#_______________RDA analysis of environmental factors and phylum————————
# 1. Loop to plot each watershed
# Load required packages
library(vegan)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(rdacca.hp) # Hierarchical partitioning required
library(gridExtra)
library(showtext)
library(grid)
# water Read data
phy <- read.delim('phy_env1252.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
sampledata <- subset(phy, type == "water", select = 11:23)
sampledata <- decostand(sampledata, method = "hellinger")
head(sampledata)
env <- subset(phy, type == "water", select = 24:33)
head(env)
env <- scale(env)
# Read group data
group <- subset(phy, select = "climate", type == "water")
# colnames(group)[colnames(group) == "NAME"] <- "sub" # Rename column
head(group)
col <- c('#E64B35', '#00A087', '#3C5488', "#FF4081", "#3498DB", "#FFC107", "#9C27B0", "#FF007F", "#87CEEB")
# '#E64B35', '#00A087', '#3C5488', "#FF4081", "#3498DB", "#FFC107", "#9C27B0"
# Perform DCA analysis (GL < 3 choose RDA, greater than 4 choose CCA, 3-4 arbitrary choice)
dca <- decorana(veg = sampledata)
dca1 <- max(dca$rproj[,1])
dca2 <- max(dca$rproj[,2])
dca3 <- max(dca$rproj[,3])
dca4 <- max(dca$rproj[,4])
GL <- data.frame(DCA1 = c(dca1), DCA2 = c(dca2), DCA3 = c(dca3), DCA4 = c(dca4))
GL
# RDA analysis
rda <- rda(sampledata, env, scale = TRUE)
r <- RsquareAdj(rda) # Raw R2
r2 <- round(r$adj.r.squared, 3) # Adjusted R2
vif.cca(rda)
rdascore <- scores(rda) # Read results
rdascore$sites # Sample scores
rda$CCA$biplot # Environmental factor response variable scores
rdascore$species # Species scores

rdascore <- scores(rda) # Read results
# Export main information
# write.csv(rdascore$sites, file = "rda_water.sample.csv")
# write.csv(rda$CCA$biplot, file = "rda_water.env.csv")
# write.csv(rdascore$species, file = "rda_water.species.csv")

# Hierarchical partitioning, analyze the contribution rate of each variable.
# Code reference https://mp.weixin.qq.com/s/L9DEzu5fc7pkzY8RtcbLlA 
library(rdacca.hp)
# res <- rdacca.hp(sampledata, env, method = "RDA", type = "R2")
# Directly use adjusted R2
res <- rdacca.hp(sampledata, env, method = "RDA", type = "adjR2")
# Save hierarchical partitioning results
write.csv(res$Hier.part, file = "water Hierarchical Partitioning.csv", row.names = TRUE)
# View results
res
# Plot
plot(res)
# View results
ph <- plot(res, lwd = 5, cex.axis = 5, main = group)
# Plot
ggsave(paste0("hp_water.jpg"), dpi = 300, plot = ph, width = 16, height = 12, units = "in")
# Permutation test, using Monte Carlo permutation test, it is possible to determine whether environmental factors are significantly related to species community,
# If the result p > 0.05, it indicates that the environmental factor is not a significant explanation for the changes in the species community,
# In other words, it is not a major influencing factor, the Monte Carlo permutation test results are as follows.
envfit <- envfit(rda, env, permutations = 999)
R <- as.matrix(envfit$vectors$r)
p <- as.matrix(envfit$vectors$pvals)
env.p <- cbind(R, p)
colnames(env.p) <- c("r2", "p-value")
KK <- as.data.frame(env.p)
KK$significance <- ifelse(KK$`p-value` < 0.001, "***", ifelse(KK$`p-value` < 0.01, "**", ifelse(KK$`p-value` < 0.05, "*", "")))
write.csv(KK, file = "water_rdaenvfit.csv") 
# Environmental factor information
RDAE <- as.data.frame(rda$CCA$biplot[, 1:2]) # Environmental factor variable scores
RDAS1 <- rdascore$sites[, 1] * 0.3 # Multiply by 0.2 to make the graph more aesthetically pleasing, does not affect the analysis
RDAS2 <- rdascore$sites[, 2] * 0.3
plotdata <- data.frame(rownames(rdascore$sites), RDAS1, RDAS2, group$climate)
colnames(plotdata) <- c("sample", "RDAS1", "RDAS2", "group")
# Add species information (multiply by 0.2 for better graph appearance)
spec <- data.frame(rdascore$species[, 1] * 0.6, rdascore$species[, 2] * 0.6)
colnames(spec) <- c("RDA1", "RDA2")
# Add label information
rda1 <- round(rda$CCA$eig[1] / sum(rda$CCA$eig) * 100, 2) # First principal axis label
rda2 <- round(rda$CCA$eig[2] / sum(rda$CCA$eig) * 100, 2) # Second principal axis label

# RDA plot绘图 
p <- ggplot(plotdata, aes(RDAS1, RDAS2)) +
  geom_point(aes(fill = group, color = group), size = 5) + 
  scale_fill_manual(values = col) +
  stat_ellipse(aes(fill = group), alpha = 0.1, color = "black", linetype = "dashed", geom = "polygon") +
  # Encircle each group with an ellipse
  xlab(paste("RDA1 (", rda1, "%)", sep = "")) + 
  ylab(paste("RDA2 (", rda2, "%)", sep = "")) +
  geom_segment(data = RDAE, aes(x = 0, y = 0, xend = RDAE[, 1], yend = RDAE[, 2]),
               colour = "red", size = 2,
               arrow = arrow(angle = 30, length = unit(0.3, "cm"))) +
  geom_text_repel(data = RDAE, segment.colour = "red",
                  aes(x = RDAE[, 1], y = RDAE[, 2], label = rownames(RDAE)), size = 10, color = "red") +
  geom_segment(data = spec
