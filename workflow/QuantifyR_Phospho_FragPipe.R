# Process ----


# Packages
library(Biostrings)
library(devtools)
library(tidyverse)
library(seqinr)
library(reshape2)


# Or install packages
#install.packages("BiocManager"); BiocManager::install("Biostrings")
#install.packages("devtools")
#install.packages("tidyverse")


# Set the working directory.
setwd("") # ???


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Process.R"))

#Load report 
report <- read_tsv("report.tsv", col_types = cols()) # Load data report


pepm <- report %>%
  remove_crap() %>% 
  select(Run, Modified.Sequence, Precursor.Charge, Stripped.Sequence, Accession = Protein.Group, `STY:79.96633`, Precursor.Normalised) %>% # Change `STY:79.96633 based on mod of interest`
  dcast(.,
        Accession + Modified.Sequence + Precursor.Charge + Stripped.Sequence + `STY:79.96633` ~ Run, 
        value.var = "Precursor.Normalised",
        fun.aggregate = mean,
        na.rm = TRUE)


samples <- 6:10 # ???


#Filter for phosphomodifications and remove other modifications from sequence
pepm2 <- pepm %>% 
  filter(str_detect(Modified.Sequence, "\\(UniMod\\:21\\)")) %>% # Change modification based on mod of interest
  mutate(Modified.Sequence = str_remove_all(Modified.Sequence, "\\(UniMod\\:1\\)")) %>% # Change modifications based on fixed and variable mods
  mutate(Modified.Sequence = str_remove_all(Modified.Sequence, "\\(UniMod\\:4\\)")) %>% 
  mutate(Modified.Sequence = str_remove_all(Modified.Sequence, "\\(UniMod\\:35\\)")) %>% 
  mutate(Feature = row_number()) 


# Load protein sequence database and parse fasta file into dataframe #
database <- read.fasta("", as.string = TRUE, seqtype = "AA", set.attributes = FALSE)  # ???
database <- database %>% 
  database_parsing(.)



pepm3 <- pepm2 %>% 
  get_identifier_FragPipe(., database, mod = "\\(UniMod\\:21\\)") %>%    # Change mod based on the modification of interest
  reduce_identifiers_FragPipe(., samples)


data <- pepm3 %>%
  select(Identifier, samples) %>%
  replace(is.na(.), 0) %>% 
  data.frame() 

# Analyze ----


# Packages
library(imp4p)
library(broom)


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))


# Workflow
# Define column indices for replicates in each condition.
a <- c(2:5)    
b <- c(6:9)    
c <- c(10:13) 
d <- c(14:17)  

# --- Group list ---
group <- list("" = a, "" = b, "" = c, "" = d) # name by condition

group.compare <- list(
  "" = list(a,b),
  "" = list(c,d),
  "" = list(b, c),
  "" = list(a, c),
  "" = list(b, d),
  "" = list(a, d)
) # name what groups to compare

# Rename the abundance columns in a simplified "Condition-Replicate" format
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
data3 <- data3 %>%
  filter(FDR < 0.05) %>%
  filter(abs(`0-60_FC`) >= 1) %>%
  calculate_hclust(., group, k = 2) %>%
  left_join(data3, ., by = names(data3)[1])


# Annotate ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Annotate.R"))


# Get Localization #
data4 <- data2 %>% 
  get_phospho_localization(., database, Threshold = 0) # Change threshold based of localization score cutoff


# Plot ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))


# PCA
# Principal component analysis (PCA)
data3 %>% plot_pca(., group) +
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
  #filter(abs(`FC_max`) >= 1) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()


# GO Summary
data4 %>%
  filter(Cluster != "NA") %>%
  plot_GO_cluster(., column = "Gene ontology", top = 3) +
  theme_custom(base_size = 24)
