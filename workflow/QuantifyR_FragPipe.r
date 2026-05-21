library(devtools)
library(dplyr)
library(tidyverse)
library(reshape2)


#set up working directory
setwd("") # ???


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Process.R"))

#Load in report.tsv file from DIANN
report <-  read_tsv("report.tsv", col_types = cols()) # Load data report


# Count total number of peptides per protein
Protm <- report %>% 
  group_by(Protein.Group) %>% 
  mutate(uniqueness = str_detect(`All Mapped Proteins`, ","), 
         `Unique Peptides` = n_distinct(Stripped.Sequence[uniqueness == FALSE])) %>% 
  mutate(Peptides = n_distinct(Stripped.Sequence)) %>% 
  ungroup()

# Remove cRAP and filter for quantifiable proteins
Protm2 <- Protm %>%
  remove_crap() %>% 
  group_by(Protein.Group) %>% 
  filter(`Peptides` >= 2 & `Unique Peptides` >= 1) %>% 
  ungroup() 

# Pivot data frmae

Protm3 <- Protm2 %>% 
  select(Run, Accession = Protein.Group, Quality = Q.Value, Peptides, `Unique Peptides`, PG.MaxLFQ) %>% 
  dcast(.,
    Accession + Peptides + `Unique Peptides` ~ Run, 
    value.var = "PG.MaxLFQ",
    fun.aggregate = mean,
    na.rm = TRUE)

Unique_Peptides <- sum(Protm3$`Unique Peptides`, na.rm = TRUE)
Peptides <- sum(Protm3$Peptides)


samples <- 1:16 #assign columns containing abundance value

data <- Protm4 %>% 
  select(Accession = Protein, samples) %>% 
  replace(is.na(.), 0) %>% 
  data.frame() 

# The rest of your QuantifyR: load in comparitive analysis tools 
library(imp4p)
library(broom)


url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))

# --- Define column indices for replicates within final data frame---
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

# --- Rename abundance columns ---
data <- data %>%
  rename_columns(., group)


# --- Clean, transform, and impute for statistical testing ---
data2 <- data %>% 
  clean_min(group, nonzero = 3) %>%
  transform_data(group, method = "log2") %>%
  impute_imp4p(group)

# --- Statistical calculations ---
# Pairwise *t*-test and one way anova test (BH FDR correction)
data3 <- data2 %>%
  calculate_ttest(., group.compare) %>% # Pairwise t-test
  calculate_1anova(., group)


# Fold change between specified comparisons
data3 <- data3 %>%
  calculate_fc(., group.compare) %>%
  add_fc_max()


# --- Clustering for abundance trends (specify number of K clusters)---
data3 <- data3 %>% 
  filter(FDR < 0.05) %>% 
  calculate_hclust(., group, k = 6) %>%
  left_join(data3, ., by = names(data3)[1])

# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))

#data4 <- data3 %>%
#add_missingness(., data, group) 

#Plotting PCA plot
data3 %>%
  plot_pca(., group) +
  theme_custom()


# Volcano Plot, this does the actual plotting
data3 %>%
  plot_volcano(.,
               group,
               group.compare,
               fdr = TRUE,
               threshold = 2,
               xlimit = 5,
               ylimit = 8) +
  theme_custom() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text( hjust = 1, size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  )
ggsave("", width = 10, height = 7, dpi = 300)

#protein abundance trends
data3 %>%
  filter(FDR < 0.05) %>% # ???
  plot_hclust(., group, k = 6) +
  theme_custom() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )

ggsave("", width = 9, height = 7, dpi = 300)



