# Setup
setwd("<YOUR_DIRECTORY_PATH>")  # Set your working directory
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
library(readxl)
library(dplyr)

# Read the data from Excel file
df <- read_excel("<PATH_TO_EXCEL_FILE>")
all_df <- df %>% filter(Group_cat != "HC")
VH <- df %>% filter(Group_cat == "VH")
VH$Group <- NULL
NOVH <- df %>% filter(Group_cat == "NOVH")
NOVH$Group <- NULL

# Everyone in one group - df remains as df

# Extract specific variables of interest from VH
VH_variables <- VH %>% select(
  FC_sig_subnetwork_FD,
  Age_at_rsfmri_scan,
  GENDER,
  LEDD,
  `PD_duration(years)`,
  TOT_UPDRS_3,
  TOT_RBD,
  alpha_syn_no_int,
  beta_amyloid_no_int,
  tau_no_int,
  Semantic_fluency_animals,
  `Hopkins Verbal Learning Test_TOTAL_RECALL`,
  `Letter-Number Sequencing_TOTRAW`,
  JLO_TOT,
  MOCATOT,
  MCA_Visuoconstr_composite,
  Moca_attention_composite
)

# Prepare dataset with added CSF measures for correlation analysis
dataset_corr <- all_df[, c(
  "Age_at_rsfmri_scan", "LEDD", "PD_duration(years)",
  "TOT_UPDRS_3", "TOT_RBD", "Semantic_fluency_animals",
  "Hopkins Verbal Learning Test_TOTAL_RECALL", "Letter-Number Sequencing_TOTRAW",
  "JLO_TOT", "MOCATOT", "alpha_syn_no_int", "beta_amyloid_no_int", "tau_no_int"
)]
col_names <- c(
  "Age", "LEDD(mg)", "PD duration (y)", "UPDRS III", "RBD Score",
  "Exec Func (sem fluency)", "Episodic Mem (HopkVerb Tot)", "Attention (Lett-num)",
  "Visuosp proc (JLO)", "MoCa Tot", "alpha_syn", "beta_amy", "tau"
)
names(dataset_corr) <- col_names

# Drop patients with missing CSF observations
dataset_corr <- na.omit(dataset_corr)

# Perform PCA
library("FactoMineR")
library("factoextra")
library("missMDA")

res.pca <- PCA(dataset_corr, graph = FALSE)

# Print PCA results
print(res.pca)

# Eigenvalues greater than 1 indicate significant variance accounted for by the PCs
eig.val <- get_eigenvalue(res.pca)
eig.val

# Plot eigenvalues
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Get more information about the PCA results
var <- get_pca_var(res.pca)
head(var$coord)  # Variable coordinates
head(var$cos2)  # Quality on the factor map (cos2)
head(var$contrib)  # Contributions to the principal components
head(var$cor)  # Correlations between variables and dimensions

# Plot variables on the factor map
fviz_pca_var(res.pca, col.var = "black")

# Plot the quality of variables on the factor map (cos2)
corrplot(var$cos2, is.corr=FALSE, tl.col="black")
# Visualize cos2 values on the factor map
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color variables by cos2 values on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, axes = c(1, 3))

# Plot correlations between variables and dimensions
corrplot(var$cor[,1:4], is.corr=FALSE, tl.col="black", col=colorRampPalette(c("#00688B","white","#EE5C42"))(300), col.lim = c(-1,1))
# Plot contributions of variables to PCs
corrplot(var$contrib[,1:4], is.corr=FALSE, tl.col="black", cl.pos='b')

# Plot contributions of variables to each of the PC 
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# Plot contributions of variables to PC1-PC4 together
fviz_contrib(res.pca, choice = "var", axes = 1:4, top = 10)

# Highlight most important variables on the correlation plot
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Extract individual values
ind <- get_pca_ind(res.pca)
# Loadings (variance-scaled or full-scaled coordinates of the variables)
loadings <- ind$coord
# Export results into TXT files
write.infile(ind, "path/to/output/ind.csv", sep = ",")
write.infile(loadings, "path/to/output/loadings.csv", sep = ",")


### Plotting correlations and heatmaps
# Load required libraries
library(ggplot2)
library(reshape2)
library(Hmisc)

# Function to format string more concisely
abbreviateSTR <- function(value, prefix) {
  lst = c()
  for (item in value) {
    if (is.nan(item) || is.na(item)) { # If item is NaN, return empty string
      lst <- c(lst, '')
      next
    }
    item <- round(item, 2) # Round to two digits
    if (item == 0) { # If rounding results in 0, clarify
      item = '<.001'
    }
    item <- as.character(item)
    item <- sub("(^[0])+", "", item)    # Remove leading 0: 0.05 -> .05
    item <- sub("(^-[0])+", "-", item)  # Remove leading -0: -0.05 -> -.05
    lst <- c(lst, paste(prefix, item, sep = ""))
  }
  return(lst)
}

# Read the data from CSV
df_all <- read.csv("~/Desktop/pca_all_loadings_CSF_72.csv")
# Select the relevant columns (2 to 5 are the PC dimensions 1-4 and 7 is the column with the average NBS FC)
d <- df_all[, c(2, 3, 4, 5, 7)]
# Compute the correlation matrix
cormatrix <- rcorr(as.matrix(d), type = 'spearman')
cormatrix$adj <- p.adjust(cormatrix$P, method = "none") ## FDR correction is later applied 

# Convert the correlation matrix to long format
cordata <- melt(cormatrix$r)
# Add labels for correlation, p-values, and adjusted p-values
cordata$labelr <- abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP <- abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj <- abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label <- paste(cordata$labelr, "\n", cordata$labelPadj, sep = "")
cordata$strike <- ""
cordata$strike[cormatrix$adj > 0.05] <- "X" # Show only significant correlations
# Set text size for plot
txtsize <- par('din')[2] / 2

# Plot the correlation heatmap for all patients
all_corr <- ggplot(cordata, aes(x = Var1, y = Var2, fill = value), type = "upper") +
  geom_tile() +
  xlab("") +
  ylab("") +
  labs(title = "PD All patients PCA loadings") +
  geom_text(label = cordata$label, size = txtsize * 1.5, color = "black") +
  geom_text(label = cordata$strike, size = txtsize * 3, color = "black", alpha = 0.5) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white",
                      midpoint = 0, limit = c(-1, 1), space = "Lab",
                      name = "Spearman\nCorrelation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1)) +
  scale_x_discrete(labels = c("PCA1", "PCA2", "PCA3", "PCA4",   "NBS FC")) +
  scale_y_discrete(labels = c("PCA1", "PCA2", "PCA3", "PCA4", "NBS FC"))
  
# Plot the correlation heatmap for VH group
dVH <- d[1:23, ]
cormatrix <- rcorr(as.matrix(dVH), type = 'spearman')
cormatrix$adj <- p.adjust(cormatrix$P, method = "none")
cordata <- melt(cormatrix$r)
cordata$labelr <- abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP <- abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj <- abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label <- paste(cordata$labelr, "\n", cordata$labelPadj, sep = "")
cordata$strike <- ""
cordata$strike[cormatrix$adj > 0.05] <- "X"
txtsize <- par('din')[2] / 2
all_corr_VH <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="PDVH PCA loadings") +
  geom_text(label=cordata$label, size=txtsize*1.5, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.5) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 10, hjust = 1)) +
  scale_x_discrete(labels = c("PCA1", "PCA2","PCA3", "PCA4","NBS FC")) + 
  scale_y_discrete(labels = c("PCA1", "PCA2", "PCA3","PCA4", "NBS FC"))
all_corr_VH

# Just looking at the pvals related to the 4 PC dimensions and NBS FC in VH from the heatmap above
pvals_VH <- c(0.02,0.67,0.37,0.39)
adjusted_p_values <- p.adjust(pvals_VH, method = "fdr") # Apply FDR correction to the p-values
# Print the original and adjusted p-values
cat("Original p-values:", pvals_VH, "\n")
cat("Adjusted p-values:", adjusted_p_values, "\n")



# NOVH
dNOVH <-d[24:72, ]
cormatrix = rcorr(as.matrix(dNOVH), type='spearman')
cormatrix$adj = p.adjust(cormatrix$P, method = "none") 
cordata = melt(cormatrix$r) 
cordata$labelr = abbreviateSTR(melt(cormatrix$r)$value, 'r')
cordata$labelP = abbreviateSTR(melt(cormatrix$P)$value, 'P non-adjusted')
cordata$labelPadj = abbreviateSTR(melt(cormatrix$adj)$value, 'P')
cordata$label = paste(cordata$labelr, "\n", 
                      cordata$labelPadj, sep = "")
cordata$strike = ""
cordata$strike[cormatrix$adj > 0.05] = "X" # changed this to only show the significant ones
txtsize <- par('din')[2] / 2
all_corr_NOVH <- ggplot(cordata, aes(x=Var1, y=Var2, fill=value), type = "upper") + geom_tile() + 
  #theme(axis.text.x = element_text(angle=90, hjust=TRUE, size=16)) +
  xlab("") + ylab("") + labs(title="PDNOVH PCA loadings") +
  geom_text(label=cordata$label, size=txtsize*1.5, color="black") + 
  geom_text(label=cordata$strike, size=txtsize*3 , color="black", alpha=0.5) +
  scale_fill_gradient2(low = "#00688B", high = "#EE5C42", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 1, 
                                   size = 10, hjust = 1)) +
  scale_x_discrete(labels = c("PCA1", "PCA2","PCA3", "PCA4","NBS FC")) + 
  scale_y_discrete(labels = c("PCA1", "PCA2", "PCA3","PCA4", "NBS FC"))
all_corr_NOVH

# Just looking at the pvals related to the 4 PC dimensions and NBS FC in NOVH from the heatmap above
pvals_NOVH <- c(0.22,0.11,0.001,0.04)
adjusted_p_values <- p.adjust(pvals_NOVH, method = "BH") # Apply FDR correction to the p-values
# Print the original and adjusted p-values
cat("Original p-values:", pvals_NOVH, "\n")
cat("Adjusted p-values:", adjusted_p_values, "\n")


## Final plots - Visualizing the relationship between NBS FC and PCA1, PCA3 and PCA4 loadings of interest
# plot correlation between NBS FC and PCA1 and separately PCA3, colouring by group
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
df_scatter <- read.csv("~/Desktop/pca_all_loadings_CSF_72.csv")
df_scatter$Group <- as.factor(df_scatter$Group)

PCA1_scatter <- ggplot(df_scatter, aes(x=Dim.1, y=NBS_FC)) + 
  geom_point(aes(alpha= 0.8,colour=Group)) +
  geom_smooth(method='lm') +
  ggtitle('NBS FC vs. PC1 loadings') +
  xlab('PC1 Loadings') + ylab('NBS FC')  + theme_minimal() +
  theme(legend.position="bottom") +  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_colour_manual(values = c("#00AFBB", "#E7B800"))
PCA1_scatter 

PCA3_scatter <- ggplot(df_scatter, aes(x=Dim.3, y=NBS_FC)) + 
  geom_point(aes(alpha= 0.8,colour=Group)) +
  geom_smooth(method='lm') +
  ggtitle('NBS FC vs. PC3 loadings') +
  xlab('PC3 Loadings') + ylab('NBS FC')  + theme_minimal() +
  theme(legend.position="bottom") +  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  scale_colour_manual(values = c("#00AFBB", "#E7B800"))
PCA3_scatter 

## Now looking at the plots separately per group separately
# Scatter plot for PDVH group
PCA1_scatter_PDVH <- ggplot(df_scatter[df_scatter$Group == "VH", ], aes(x = Dim.1, y = NBS_FC)) + 
  geom_point(aes(alpha = 0.8, colour = Group)) +
  geom_smooth(method = 'lm') +
  ggtitle('NBS FC vs. PC1 loadings - PDVH') +
  xlab('PC1 Loadings') + ylab('NBS FC') + theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#00AFBB")) +
  scale_colour_manual(values = c("#00AFBB"))
PCA1_scatter_PDVH

PCA3_scatter_PDVH <- ggplot(df_scatter[df_scatter$Group == "VH", ], aes(x = Dim.3, y = NBS_FC)) + 
  geom_point(aes(alpha = 0.8, colour = Group)) +
  geom_smooth(method = 'lm') +
  ggtitle('NBS FC vs. PC3 loadings - PDVH') +
  xlab('PC3 Loadings') + ylab('NBS FC') + theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#00AFBB")) +
  scale_colour_manual(values = c("#00AFBB"))
PCA3_scatter_PDVH

# Scatter plot for PDNOVH group
PCA1_scatter_PDNOVH <- ggplot(df_scatter[df_scatter$Group == "NOVH", ], aes(x = Dim.1, y = NBS_FC)) + 
  geom_point(aes(alpha = 0.8, colour = Group)) +
  geom_smooth(method = 'lm') +
  ggtitle('NBS FC vs. PC1 loadings - PDNOVH') +
  xlab('PC1 Loadings') + ylab('NBS FC') + theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#E7B800")) +
  scale_colour_manual(values = c("#E7B800"))
PCA1_scatter_PDNOVH

PCA3_scatter_PDNOVH <- ggplot(df_scatter[df_scatter$Group == "NOVH", ], aes(x = Dim.3, y = NBS_FC)) + 
  geom_point(aes(alpha = 0.8, colour = Group)) +
  geom_smooth(method = 'lm') +
  ggtitle('NBS FC vs. PC3 loadings - PDNOVH') +
  xlab('PC3 Loadings') + ylab('NBS FC') + theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#E7B800")) +
  scale_colour_manual(values = c("#E7B800"))
PCA3_scatter_PDNOVH

PCA4_scatter_PDNOVH <- ggplot(df_scatter[df_scatter$Group == "NOVH", ], aes(x = Dim.4, y = NBS_FC)) + 
  geom_point(aes(alpha = 0.8, colour = Group)) +
  geom_smooth(method = 'lm') +
  ggtitle('NBS FC vs. PC4 loadings - PDNOVH') +
  xlab('PC4 Loadings') + ylab('NBS FC') + theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#E7B800")) +
  scale_colour_manual(values = c("#E7B800"))
PCA4_scatter_PDNOVH
