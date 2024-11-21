setwd("E:/Rpractice/meta/tem_article/randomforest_meta/water")
#R language | Calculate microbial community diversity and its environmental drivers
library(ggplot2)  
library(vegan)
library(randomForest)
library(rfPermute)
####Correlation analysis between environmental variables and species abundance
library(psych)
library(reshape2)
library(patchwork)#Merge images
#Link: https://pan.baidu.com/s/1SSSKA6eB-6Efex2J4BumGQ  Access code: isqe
#Data according to the corresponding OTU table, species file, and grouping file
phy<-read.delim('phy_env_water903_alt3_new.txt', sep = '\t', row.names = 1, check.names = FALSE)

df<-subset(phy,type=="water",c("shannon_diversity","WT","DO"  ,"pH","TN"  ,"TP","NH4_N","PRE","TEM","alt_rst"))

# Set a random seed to make the results reproducible
set.seed(1234)
# Use the randomForest package to calculate the contribution importance of environmental factors
rf_results <- rfPermute(shannon_diversity~ ., data = df[, -2], importance = TRUE, ntree = 1000, num.cores = 15)
# Extract the explanatory power of the predictive factors
predictor_var <- data.frame(importance(rf_results, scale = TRUE), check.names = FALSE)
# Extract the significance of the predictive variables
predictor_sig <- as.data.frame((rf_results$pval)[,,2])
# Set column names
colnames(predictor_var) <- c("var", "p", "IncNodePurity", "IncNodePurity.p")
# Set flags based on significance
df_pre <- predictor_var
df_pre$sig[df_pre$p < 0.05] <- "*"
df_pre$sig[df_pre$p >= 0.05] <- " "
# Set the color for significant factors
df_pre$IncNodePurity[df_pre$sig == "*"] <- "orange"
# Set the color for non-significant factors
df_pre$IncNodePurity[df_pre$sig == " "] <- "skyblue"
# Sort by importance
#df_pre <- df_pre[order(df_pre$var, decreasing = TRUE), ]
# Save the data
write.csv(df_pre,file=paste0("water environmental factors explaining diversity.csv"))
# Plot a bar chart
# Custom order list
custom_order <- c("WT","DO"  ,"pH","TN"  ,"TP","NH4_N","PRE","TEM",	"alt_rst")
# Create a row name column
df_pre$rowname <- factor(rownames(df_pre), levels = custom_order)
df_pre <- df_pre[order(df_pre$rowname), ]#Sort by specified order
p1 <- ggplot(data = df_pre, aes(y= rowname, x= df_pre$var)) +
  geom_bar(stat = 'identity', width = 0.8, fill = df_pre$IncNodePurity) +
  theme_bw() +#coord_flip()+#Flip x and y axes
  labs(x = 'Increase in MSE (%)', y = '') +scale_x_continuous(expand = c(0,0)) +
  geom_text(aes(label = df_pre$sig, x = df_pre$var + 0.28, y = rowname), vjust =-0.35, size = 9) +
  theme_bw(base_line_size = 1.05, base_rect_size = 1.05) +theme(axis.text.y = element_blank())+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9, size = 11), 
        axis.title.x.top = element_text(vjust = 0.5)) +
  ggtitle("water")
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
#theme(axis.text = element_text(color = "black", size = 11))

#Read the environmental variable and species abundance matrix
env <- df$shannon_diversity
spe <- df[,-(1:2)]
#spe <- spe[rownames(env), ]

#Can be executed via the psych package function corr.test()
#Here, Pearson's correlation coefficient is used, and no p-value correction is made for now (can be specified via the adjust parameter)
pearson <- corr.test(env, spe, method = 'pearson', adjust = 'none')
r <- data.frame(pearson$r)  #Pearson correlation coefficient matrix
p <- data.frame(pearson$p)  #P-value matrix
#Result organization for plotting
r$env <- rownames(r)
p$env <- rownames(p)
r <- melt(r, id = 'env')
p <- melt(p, id = 'env')
pearson <- cbind(r, p$value)
colnames(pearson) <- c('env', 'spe', 'pearson_correlation', 'p.value')
pearson$spe <- factor(pearson$spe, levels = colnames(spe))
head(pearson)  #Organized table of Pearson correlation statistics between environmental variables and species abundance

#Use ggplot2 to plot a heatmap of Pearson correlation between environmental variables and species abundance
# Custom order list
custom_order <- c("WT","DO"  ,"pH","TN"  ,"TP","NH4_N","PRE","TEM",	"alt_rst")

#Convert 'env' column to factor and sort according to custom order
pearson$env <- factor(pearson$env, levels = custom_order)
#If you want to mark the significance of the Pearson correlation coefficient on the plot, refer to the following operations
pearson[which(pearson$p.value<0.001),'sig'] <- '***'
pearson[which(pearson$p.value<0.01 & pearson$p.value>0.001),'sig'] <- '**'
pearson[which(pearson$p.value<0.05 & pearson$p.value>0.01),'sig'] <- '*'
head(pearson)  #Organized table of Pearson correlation statistics between environmental variables and species abundance
#Plot the heatmap using the sorted data
p2 <- ggplot() +
  geom_tile(data = pearson, aes(x = env, y = spe, fill = pearson_correlation)) +
  scale_fill_gradientn(colors = c('#00B554', 'white', '#B586CA'), limit = c(-1, 1)) +
  theme(panel.grid = element_line(), panel.background = element_rect(color = 'black'), 
        legend.key = element_blank(), legend.position = "bottom",
        #legend.margin = margin(t = -1, unit = "cm"),  # Adjust the spacing between the legend and the top of the plot
        #legend.box.margin = margin(t = 0, unit = "cm"),  # Adjust the spacing between the legend contents and the top
        axis.text.x = element_text(color = 'black', angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio=1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(y = '', x = '', fill = '')+
  geom_text(data = pearson, aes(x = env, y = spe, label = sig), size = 8)

p3<-p2+p1 
# Add title

ggsave(paste0("water environmental factors explaining diversity significance.pdf"),dpi = 300,plot = p3, width = 8, height = 12, units = "in")


#——————————————————Loop____
data<-subset(phy,type=="water",c("shannon_diversity","alt3","WT","DO"  ,"pH","TN"  ,"TP","NH4_N","PRE","TEM",	"alt_rst"))
group<-unique(data$climate)
library(ggplot2)
library(gridExtra)
groups<-c("S","SW","C", "N","QT")
groups<-c("L","M"