# Set working directory and load necessary libraries
setwd("E:/Rpractice/meta/tem_article/distance_climate")
library(vegan)
library(ggplot2)
library(reshape2)
library(geosphere)
library(gridExtra)

# Read data
phy <- read.delim('phy_env_water903_alt3_new.txt', sep = '\t', row.names = 1, check.names = FALSE)
water_phylum <- subset(phy, type == "water")
site <- water_phylum[, c("longitude", "latitude")]
site_dis <- geosphere::distm(site[c("longitude", "latitude")])
site_dis <- site_dis / 1000  # Convert distance to km
rownames(site_dis) <- rownames(site)
colnames(site_dis) <- rownames(site)
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)

# Extract environmental data
water <- water_phylum[, c(11:23)]
comm_sim <- 1 - as.matrix(vegan::vegdist(water, method = 'bray'))
diag(comm_sim) <- 0
comm_sim[upper.tri(comm_sim)] <- 0
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)

# Merge distance and community similarity data
comm_dis <- merge(comm_sim, site_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'site_dis')
write.table(comm_dis, 'water_distance_decay.csv', sep = '\t', row.names = FALSE, quote = FALSE)

# Plotting
theme_set(theme_bw())
a <- ggplot(comm_dis, aes(x = site_dis, y = comm_sim)) + 
  geom_point(size = 0.15, color = "skyblue") +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, color = "grey", fullrange = TRUE) +
  xlab("Geographical Distance (km)") + 
  ylab("Bray-Curtis Similarity") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.line = element_line(color = "black", size = 0.75), 
        axis.title = element_text(size = 16), 
        panel.border = element_rect(color = "black", size = 0.75)) +
  ggtitle("Water")

# Linear regression
fit <- lm(comm_sim ~ site_dis, data = comm_dis)
summary(fit)

# Get p-value and slope
p_value <- summary(fit)$coefficients[2, 4]
slope <- fit$coefficients[2]

# Construct p-value label
p_label <- ifelse(p_value < 0.001, "p < 0.001", paste("p =", formatC(p_value, digits = 3, format = "f")))
spearman_corr <- cor(comm_dis$site_dis, comm_dis$comm_sim, method = "spearman")

# Add R-squared and p-value text
a <- a +
  annotate("text", x = Inf, y = Inf, label = paste("R2 =", round(summary(fit)$r.squared, 4)), hjust = 1, vjust = 1, size = 6) +
  annotate("text", x = Inf, y = Inf, label = p_label, hjust = 1, vjust = 2, size = 6) +
  annotate("text", x = Inf, y = -Inf, label = paste("N =", nrow(comm_dis)), hjust = 1, vjust = 0, size = 5) +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", round(spearman_corr, 3)), hjust = 1, vjust = -1, size = 6)

# Save plot
ggsave("water_distance_decay.pdf", plot = a, dpi = 300, width = 4, height = 4, units = "in")

# Repeat the above steps for sediment data and combine plots