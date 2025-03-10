# This code plots the NBS Subnetwork FC in ICICLE, which is the masked version of PPMI.

# Load the necessary libraries.
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
library(readxl)
library(dplyr)
library(ggpubr)
library(MetBrewer)
library(viridis)
library(ggsci) # includes nature journal palettes
library(jcolors)
library(tidyr)

# Set the working directory.
setwd("/path/to/directory")

# Read the data.
df <- read_excel("dataset.xlsx", sheet = "data")

# Convert the `Group` column to a factor.
df$Group <- as.factor(df$Group)

# Filter the data to only include participants from the ICICLE study.
df <- df %>% filter(Study == "ICICLE")

# Create a new column called `ICICLE_NBS_masked_FC`.
df$ICICLE_NBS_masked_FC <- df$ICICLE_NBS_masked_FC_from_PPMI

# Create a new column called `NBS_full_67edges_masked`.
df$NBS_full_67edges_masked <- df$NBS_full_67edges_masked

# Create a new column called `UPDATED_severity`.
df$UPDATED_severity <- df$UPDATED_severity

# Fit a linear model to the data.
linearModel <- lm(UPDATED_severity ~ ICICLE_NBS_masked_FC, data = df)
# Summary of the linear model.
summary(linearModel)

# Plot the data with a regression line.
ggplot(df, aes(x = ICICLE_NBS_masked_FC, y = UPDATED_severity)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ylab("Hallucination Severity") +
  xlab("Masked FC from PPMI NBS") 

# Plot the data with a regression line.
ggplot(df, aes(y = TOT_frequency, x = ICICLE_NBS_masked_FC)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ylab("Hallucination Frequency score") +
  xlab("Masked FC from PPMI NBS") 
