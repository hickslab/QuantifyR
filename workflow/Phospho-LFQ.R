# PhosphoLFQ
# Phosphosite-level label-free quantification workflow


# Process ----


# Packages
library(Biostrings)
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("BiocManager"); BiocManager::install("Biostrings")
#install.packages("devtools")
#install.packages("tidyverse")


# Set the working directory.
# The example below is based on downloading the package from GitHub.
setwd("~/ProgenesisLFQ/data") # ???


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/Process.R"))


# Load Peptide Measurements spreadsheet exported from Progenesis.
pepm <- "20190715_EWM_TOR1_Phospho_pepm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols())


# Load Peptide Measurements spreadsheet exported from Progenesis.
protm <- "20190715_EWM_TOR1_Phospho_protm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols()) %>%
  select(1:4) %>%
  separate_rows(., Accession, sep = ";")


# Load protein sequence database.
database <- "Cr_uniprot_crap_20190130.fasta" %>% # ???
  readAAStringSet(.)

names(database) <- str_split(names(database), " ", simplify = TRUE)[, 1]


# Define normalized abundance columns in pepm
samples <- 19:38 # ???


# Process data
pepm2 <- pepm %>%
  filter(Score > 13) %>%
  filter(Description != "cRAP") %>%
  left_join(., protm, by = "Accession") %>%
  reduce_features() %>%
  filter(str_detect(Modifications, "Phospho")) %>%
  get_identifier(., database, mod = "Phospho") %>%
  reduce_identifiers(., samples)

data <- pepm2 %>%
  select(Identifier, samples) %>%
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
a <- 2:6; b <- 7:11; c <- 12:16; d <- 17:21 # ???

group <- list("Control" = a, "AZD8055" = b, "Torin1" = c, "Rapamycin" = d) # ???

group.compare <- list("Control-AZD8055" = list(a, b),
                      "Control-Torin1" = list(a, c),
                      "Control-Rapamycin" = list(a, d)) # ???


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
  calculate_ttest(., group.compare, fdr = FALSE)


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
               fdr = FALSE,
               threshold = 2,
               xlimit = 8,
               ylimit = 5) +
  theme_custom()


# Trend Profiles
data4 %>%
  filter(FDR < 0.05) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()


