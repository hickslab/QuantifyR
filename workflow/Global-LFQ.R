# ProteinLFQ
# Protein-Level Label-Free Quantification


# Process ----


# Packages
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")


# Set the working directory.
# The example below is based on downloading the package from GitHub.
setwd("~/ProgenesisLFQ/data") # ???


# Load Peptide Measurements spreadsheet exported from Progenesis.
protm <- "20180502_WOS52_Cr_UPS_protm.csv" %>%
  read_csv(., skip = 2, col_types = cols())


# Separate leading protein accession from group members.
protm <- protm %>%
  separate(Accession,
           into = c("Accession", "Group members"),
           sep = ";",
           extra = "merge",
           fill = "right")


# Define normalized abundance columns in pepm
samples <- 11:22 # ???


# Process data.
protm2 <- protm %>%
  filter(Description != "cRAP") %>%
  filter(`Peptide count` >= 2 & `Unique peptides` >= 1)

data <- protm2 %>%
  select(Accession, samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_Process.csv") # ???


# Analyze ----


# Packages
library(imp4p)
library(broom)


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/Analyze.R"))


# Workflow
# Define the column indeces for replicates in each condition.
a <- 2:5; b <- 6:9; c <- 10:13 # ???

group <- list("25" = a, "50" = b, "100" = c) # ???

group.compare <- list("25-50" = list(a, b),
                      "25-100" = list(a, c),
                      "50-100" = list(b, c)) # ???


# Rename the abundance columns in a simplified "Condition-Replicate" format.
data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns
data2 <- data %>%
  clean_min(., group, nonzero = 3) %>% # ???
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)


# Pairwise *t*-test
data3 <- data2 %>%
  calculate_ttest(., group.compare, fdr = TRUE)


# One-way ANOVA
data4 <- data2 %>%
  calculate_1anova()


# Fold change
data3 <- data3 %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_fc_max()

data4 <- data4 %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_fc_max()


# Annotate ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/Annotate.R"))


# UniProt
data5 <- data4 %>%
  add_missingness(., data, group)


# Plot ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/Plot.R"))


# PCA
# Principal component analysis (PCA)
data2 %>% plot_pca(., group) +
  theme_custom()


# Volcano Plot
data3 %>%
  plot_volcano(.,
               group,
               group.compare,
               fdr = TRUE,
               threshold = 2,
               xlimit = 8,
               ylimit = 5) +
  theme_custom()


# Trend Profiles
data4 %>%
  filter(FDR < 0.05) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()
