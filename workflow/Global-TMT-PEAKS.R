# Protein
# Protein-Level TMT Quant


# Process ----


# Packages
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")


# Set the working directory.
# The example below is based on downloading the package from GitHub.
setwd("setworkingdirectory") # ???


# Load Protein Measurements spreadsheet exported from PEAKS.
protm <- "proteins_TMT_real.csv" %>% # ???
  read_csv(., col_types = cols())


# Separate leading protein accession from group members.
protm <- protm %>%
  separate(Accession,
           into = c("Accession", "Group members"),
           sep = ";",
           extra = "merge",
           fill = "right")


# Define normalized abundance columns in pepm
samples <- 7:18 # ???


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
a <- 2:4; b <- 5:7; c <- 8:10; d <- 11:13



group <- list("128" = a, "129" = b, "130" = c, "131" =d) # ???

group.compare <- list("128-129" = list(a, b),
                      "128-130" = list(a, c),
                      "128-131" = list(a, d),
                      "129-130" = list(b, c),
                      "129-131" = list(b, d),
                      "130-131" = list(c, d))# ???


# Rename the abundance columns in a simplified "condition-replicate" format.
data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns
data2 <- data %>%
  clean_min(., group, nonzero = 2) %>%  
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
