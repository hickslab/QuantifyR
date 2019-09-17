# RedoxLFQ
# Reversibly Oxidized Cys-Level Label-Free Quantification


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
setwd("~/QuantifyR/data") # ???


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Process.R"))


# Load Peptide Measurements spreadsheet exported from Progenesis.
pepm <- "20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_pepm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols())


# Load Peptide Measurements spreadsheet exported from Progenesis.
protm <- "20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_protm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols()) %>%
  select(1:4) %>%
  separate_rows(., Accession, sep = ";")


# Load protein sequence database.
database <- "Cr_uniprot_crap_20190130.fasta" %>% # ???
  readAAStringSet(.)

names(database) <- str_split(names(database), " ", simplify = TRUE)[, 1]


# Define normalized abundance columns in pepm
samples <- 19:34 # ???


# Process input data.
pepm2 <- pepm %>% 
  filter_score(., 13) %>%
  filter(Description != "cRAP") %>%
  left_join(., protm, by = "Accession") %>%
  reduce_features() %>%
  filter_redox(., reduced = "Nethylmaleimide") %>%
  get_identifier_redox(., database, reduced = "Nethylmaleimide") %>%
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
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))


# Workflow
# Define the column indeces for replicates in each condition.
a <- 2:5; b <- 6:9; c <- 10:13; d <- 14:17

group <- list("0" = a, "15" = b, "30" = c, "60" = d) # ???

group.compare <- list("0-15" = list(a, b),
                      "0-30" = list(a, c),
                      "0-60" = list(a, d)) # ???


# Rename the abundance columns in a simplified "Condition-Replicate" format.
data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns
data2 <- data %>%
  clean_min(., group, nonzero = 3) %>% # ???
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)


# Hypothesis testing
data3 <- data2 %>%
  calculate_ttest(., group.compare) %>% # Pairwise t-test
  calculate_1anova(., group) # One-way ANOVA


# Fold change
data3 <- data3 %>%
  calculate_fc(., group.compare) %>%
  add_fc_max()


# Heirarchical clustering
data3 <- data3 %>%
  filter(FDR < 0.05) %>%
  filter(abs(`0-60_FC`) >= 1) %>%
  calculate_hclust(., group, k = 2) %>%
  left_join(data3, ., by = names(data3)[1])


# Annotate ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Annotate.R"))


# Annotations
uniprot <- "Cr_uniprot_20190130_annotation.tsv"


# Only missing values
data4 <- data3 %>%
  add_missingness(., data, group)


# Accessions and annotate
data4 <- data4 %>%
  add_accession() %>%
  add_uniprot(., uniprot)


# Plot ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))


# PCA
# Principal component analysis (PCA)
data4 %>% plot_pca(., group) +
  theme_custom()


# Volcano Plot
data4 %>%
  plot_volcano(.,
               group,
               group.compare,
               fdr = TRUE,
               threshold = 2,
               xlimit = 8,
               ylimit = 5) +
  theme_custom()


# Trend Profiles
data4[group %>% flatten_int()] <- 2^data4[group %>% flatten_int()]

data4 %>%
  filter(FDR < 0.05) %>% # ???
  filter(abs(`0-60_FC`) >= 1) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()


# GO Summary
data4 %>%
  filter(Cluster != "NA") %>%
  plot_GO_cluster(., column = "Gene ontology", top = 3) +
  theme_custom(base_size = 24)

  