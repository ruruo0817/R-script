# Generalized Additive Models
# Install the mgcv package if not already installed
# install.packages("mgcv")
# Load the mgcv package for fitting Generalized Additive Models (GAMs)
library(mgcv)

# Read species abundance and water body environmental data
setwd("E:/Rpractice/meta/tem_article/GAM_alt")
phy <- read.delim('phy_env_water_alt3_new.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
# The following line seems to read the same file again, which might be an error
# phy <- read.delim('water1573.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
df <- subset(phy, type == "water", c("shannon_diversity", "alt_rst"))

# Fit a GAM model with a smooth term for alt_rst
gam_model <- gam(shannon_diversity ~ s(alt_rst, k = 5), data = df, family = gaussian()) # k=5 controls the flexibility of the fit
summary(gam_model)
gam_model
plot(gam_model, pch = 20, shade = TRUE, residuals = TRUE)
# Save the plot as a PDF or other formats
pdf("gam_model_Water2.pdf", width = 4, height = 4)
plot(gam_model, pch = 20, shade = TRUE, residuals = TRUE)
title("Water")
xlabel <- expression(paste("Altitude (", m, ")"))
ylabel <- "Shannon Diversity"
# Restore default graphical parameters
par(lwd = 1)  # Restore default line width
dev.off()

# Correlation between different bacterial phyla and altitude
phy <- read.delim('water1573.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(phy)
df <- subset(phy, type == "water", c(11:23, 33))

head(df)
# Handle outliers if there are negative or non-integer values
df <- df[df$alt_rst >= 0, ]  # Assuming alt_rst column should not be negative
# Load necessary packages
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming df is the data frame you provided
# df contains the abundance of each bacterial phylum and altitude data

# Convert data to long format
df_long <- df %>%
  gather(key = bacteria, value = abundance, -alt_rst)
# Define a function to fit GAM models and predict values
fit_and_predict_gam <- function(data){
  gam_model <- gam(abundance ~ s(alt_rst, k = 3), data = data, family = gaussian()) # k=3 controls the flexibility of the fit
  pred_data <- data.frame(alt_rst = seq(min(data$alt_rst), max(data$alt_rst), length.out = 100))
  predictions <- predict(gam_model, newdata = pred_data, type = "response")
  pred_data$abundance <- predictions
  return(pred_data)
}

# Use the function to fit and predict for long format data
pred_data <- df_long %>%
  group_by(bacteria) %>%
  do(fit_and_predict_gam(.))

# Get summaries for each model
unique_bacteria <- unique(pred_data$bacteria)

model_summaries <- lapply(unique_bacteria, function(bac) {
  bac_data <- df_long %>% filter(bacteria == bac)
  model <- gam(abundance ~ s(alt_rst, k = 3), data = bac_data, family = gaussian())
  summary(model)
})
# View results
model_summaries
# Create an empty data frame
results_df <- data.frame()

# Loop through each model's results
for (summary in model_summaries) {
  # Convert each model's results to a data frame and add to results_df
  summary_df <- data.frame(t(summary))
  results_df <- rbind(results_df, summary_df)
}

# Print results
print(results_df)

# Create a ggplot2 chart
colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
  "#aec7e8", "#ffbb78", "#98df8a"
)

p1 <- ggplot(pred_data, aes(x = alt_rst, y = abundance, color = bacteria)) +
  geom_line() +
  labs(title = "Water",
       x = "Elevation(m)",
       y = "Relative Abundance(%)") +
  scale_color_manual(values = colors) + # Use the Dark2 color palette common in Nature
  theme_classic() + # White background
  theme(
    axis.line = element_line(size = 0.5, color = "black"),
    axis.title = element_text(size = 7.5),
    axis.text = element_text(size = 7.5),
    plot.title = element_text(size = 7.5, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border  # Place legend on the right side
    legend.key.size = unit(0.1, "in"),  # Adjust the size of the legend keys
    legend.text = element_text(size = 6),  # Adjust the size of the legend text
    legend.title = element_text(size = 8, face = "bold"),
    axis.ticks = element_line(size = 0.5)  # Add axis tick marks # Adjust the size of the legend title
    )+
    expand_limits(x = 0, y = 0) +  # Set the axis to start from zero
  scale_x_continuous(expand = c(0, 0)) +  # Do not expand the axis range
  scale_y_continuous(expand = c(0, 0))  # Do not expand the axis range
p1

ggsave("WATER_GAM.tif", plot = p1, dpi = 300, width = 6, height = 5, units = "in")

#————————————Sediment————————————————
phy <- read.delim('sediment.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

df <- subset(phy, type == "sediment", c("shannon_diversity", "alt_rst"))

gam_model1 <- gam(shannon_diversity ~ s(alt_rst, k = 5), data = df, family = gaussian()) # k=5 controls the flexibility of the fit
summary(gam_model1)
gam_model1
plot(gam_model1, pch = 20, shade = TRUE, residuals = TRUE)
# Save the plot as a PDF or other formats
pdf("gam_model_sediment2.pdf", width = 4, height = 4)
plot(gam_model1, pch = 20, shade = TRUE, residuals = TRUE)
title("Sediment")
xlabel <- expression(paste("Altitude (", m, ")"))
ylabel <- "Shannon Diversity"
# Restore default graphical parameters
par(lwd = 1)  # Restore default line width
dev.off()

# Loop
phy <- read.delim('sediment1063.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

head(phy)

df <- subset(phy, type == "sediment", c(11:23, 33))
df <- subset(phy, type == "sediment", c("Proteobacteria", "Actinobacteriota", "Bacteroidota", 
                                        "Acidobacteriota", "Chloroflexi", "Firmicutes", "Cyanobacteria", "alt_rst"))

head(df)
# Handle outliers if there are negative or non-integer values
df <- df[df$alt_rst >= 0, ]  # Assuming alt_rst column should not be negative
# Load necessary packages
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
# Assuming df is the data frame you provided
# df contains the abundance of each bacterial phylum and altitude data

# Convert data to long format
df_long <- df %>%
  gather(key = bacteria, value = abundance, -alt_rst)

# Define a function to fit GAM models and predict values
fit_and_predict_gam <- function(data){
  gam_model <- gam(abundance ~ s(alt_rst, k = 3), data