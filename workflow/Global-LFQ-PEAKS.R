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
setwd("Z:/Lab_Members/Sam/Projects/CCAMP1/20220501_ztof/LFQ_script") # ???


# Load Protein Measurements spreadsheet exported from PEAKS.
protm <- "proteins.csv" %>% # ???
  read_csv(., col_types = cols())


# Separate leading protein accession from group members.
protm <- protm %>%
  separate(Accession,
           into = c("Accession", "Group members"),
           sep = ";",
           extra = "merge",
           fill = "right")


# Define normalized abundance columns in pepm
samples <- 10:29 # ???


# Process input data.
protm2 <- protm %>%
  filter(Description != "cRAP") %>%
  filter(`#Peptides` >= 2 & `#Unique` >= 1)

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
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))



# Workflow
# Define the column indices for replicates in each condition.

######### In this step, you should define the column indicies from the new data matrix, 
######### rather than the original imported csv file. 


##### Change to this #####
a <- 2:6; b <- 7:11; c <- 12:16; d <- 17:21



group <- list("WT" = a, "R1" = b, "R2" = c, "R3" =d) # ???

group.compare <- list("WT-R1" = list(a, b),
                      "WT-R2" = list(a, c),
                      "WT-R3" = list(a, d)) # ???


# Rename the abundance columns in a simplified "condition-replicate" format.
data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns
data2 <- data %>%
  clean_min(., group, nonzero = 3) %>%  
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

##### this function gives an error if you run it if there are not > 2 proteins with FDR values < 0.05
data3 <- data3 %>%
  filter(FDR < 0.05) %>%
  calculate_hclust(., group, k = 2) %>%
  left_join(data3, ., by = names(data3)[1])


# Annotate ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Annotate.R"))


# UniProt
data4 <- data3 %>%
  add_missingness(., data, group)


# Plot ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))


# PCA
data4 %>%
  plot_pca(., group) +
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
data4 %>%
  filter(FDR < 0.05) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()
