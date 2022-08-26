## TMT Peptide level analysis ##

## Lines of code marked by `# ???` at the end indicate places where the user 
## should adjust code to fit their dataset

library(devtools)
library(tidyverse)
library(Biostrings)
library(seqinr)

#install.packages("seqinr")

# set working directory #

setwd("setworkingdirectory") # ???


# Functions #

url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Process.R"))

# read in TMT protein-supporting peptides from PEAKS #

pept <- "protein-peptides.csv" %>%  # ???
  read_csv(., col_types = cols()) %>%
  select(., 1:2, "Accession" = `Protein Accession`, 4:42) # ??? change 42 to the total col length of the dataframe


# Load TMT quantified protein measurements from PEAKS #

protm <- "proteins.csv" %>%  # ???
  read_csv(., col_types = cols())

# Load fasta file and convert to data frame #

p1 <- read.fasta("inputfastafile.fasta", as.string = TRUE, seqtype = "AA", set.attributes = FALSE)  # ???

as.data.frame(p1)-> df

df2 <- df %>% 
  t() %>% 
  as.data.frame()

df2 <- df2 %>%
 rownames_to_column(var = "Accession") %>% 
 separate(Accession, 
          into = c("sp", "Accession"),
          by = "tr.") %>%  # ??? change the separation by = "" based on the leading characters before the database accessions
 select(1:2, "ProteinSequence" = V1)


# Joining proteins and filter identical peptides with different accessions #

# Filter peptides associated with the same feature #

pept2 <- pept %>% 
  left_join(., protm, by = "Accession") %>%
  filter(Description != "cRAP") %>%
  group_by(scan, Peptide, PTM.x,`-10lgP.x`) %>%
    top_n(1, `#Unique`) %>%
    top_n(1, `-10lgP.y`) %>%
    dplyr::slice(1) %>%
  
  # Summarize identical features with multiple peptides #
  group_by(scan) %>%
  top_n(1, `-10lgP.x`) %>%
  dplyr::slice(1) %>%
  ungroup()


# Filter out modification masses on the peptide sequence #

pept3 <- pept2 %>% 
  mutate(Sequence = str_remove_all(Peptide, "(\\()(.*?)(\\))"))

# Filter out residue located before and after sequence #

pept3 <- pept3 %>% 
  mutate(Sequence = str_remove(Sequence, "([K][.]|[R][.]|[M][.])?")) %>% 
  mutate(Sequence = str_remove(Sequence, "([.][A-Z])"))
  

# Filter for d0 and d5 modifications and determine Cys identifier(s) #

pept4 <- pept3 %>% 
  separate_rows(AScore, sep = "\\;") %>% 
  filter(str_detect(AScore, "N-ethylmaleimide")) %>% 
  mutate(position = str_extract(AScore, "(?<=\\C)(.*)(?=\\:[A-Z])") %>% 
           as.numeric()) %>% 
  mutate(residue = str_sub(Sequence, start = position, end = position))


# Get modification site and build identifier name #

pept5 <- pept4 %>% 
  left_join(., df2, by = "Accession") %>% 
  mutate(locate = str_locate(ProteinSequence, Sequence)[,1]) %>% 
  mutate(Identifier = (locate + position - 1) %>% 
           paste("C", ., sep = "") %>% 
           #paste(., collapse = "-")
           paste(Accession, ., sep = "--"))
      
# Build NEM modification type into identifier #

pept6 <- pept5 %>% 
  mutate( `Identifier1` = if_else((str_detect(AScore, "D5")), 
                                  paste(Identifier, "H", sep = "-"),
                                  paste(Identifier, "L", sep = "-")))


# Filter by single modification or two modifications on any one Cys identifier #

# Filter by modification type #

pept7 <- pept6 %>% 
  mutate(`HorL` = if_else((str_detect(AScore, "D5")),
                                2, ## 2 = heavy NEM modification
                                1)) %>% ## 1 = light NEM modification
  mutate(Modification = str_remove(AScore, "([:][0-9]*)([.][0-9]*)")) %>% 
  group_by(Identifier) %>% 
  mutate(mean_modification = mean(`HorL`)) %>% 
  ungroup()


# Identify those that have two or more modifications on the same Cys identifier #

pept8 <- pept7 %>% 
  mutate(`RedoxState` = if_else(1 < mean_modification & mean_modification < 2 | mean_modification > 2,
                                2,
                                1))

# Define normalized peptide abundance columns in pept2 #

samples <- 17:22 # ??? - select the normalized abundance values from the pept2 dataframe

# Reduce identifiers to only include one identifier for any one Cys #

pept9 <- pept8 %>% 
  group_by(Identifier1) %>%
  top_n(1,`-10lgP.x` ) %>%
  dplyr::slice(1) %>%
  ungroup()

pept10 <- pept8 %>% 
  select(Identifier1, samples) %>% 
  gather(sample, abundance, -1) %>% 
  mutate(sample = factor(sample, levels = pept8[, samples] %>% names())) %>% 
  
  group_by(Identifier1, sample) %>%
  summarize(sum = sum(abundance)) %>% 
  
  spread(sample, sum) %>%
  ungroup()

pept9[, samples] <- pept10[, -1]


# Find number of reduced Cys peptides #

Reduced <- pept9 %>% 
  filter(HorL == 1)

# Find number of oxidized Cys peptides #

Oxidized <- pept9 %>% 
  filter(HorL == 2)


data <- pept9 %>%
  select(Identifier1, samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_Process.csv") # ???


# Analyze Cys identifiers with only one redox state ----


# Packages #

library(imp4p)
library(broom)


# Functions #
 
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))


# Workflow #
# Define the column indeces for replicates in each condition #

a <- 2:4; b <- 5:7

group <- list("WT" = a, "H2O2" = b) # ???

group.compare <- list("WT-H2O2" = list(a, b)) # ???


# Rename the abundance columns in a simplified "Condition-Replicate" format #

data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns #

data2 <- data %>%
  clean_min(., group, nonzero = 2) %>% # ???
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)


# Hypothesis testing #

data3 <- data2 %>%
  calculate_ttest(., group.compare) %>% # Pairwise t-test
  calculate_1anova(., group) # One-way ANOVA


# Fold change #

data3 <- data3 %>%
  calculate_fc(., group.compare) %>%
  add_fc_max()



# Heirarchical clustering #

data3 <- data3 %>%
  filter(FDR < 0.05) %>%
  filter(abs(`0-60_FC`) >= 1) %>%
  calculate_hclust(., group, k = 2) %>%
  left_join(data3, ., by = names(data3)[1])

# Plot Cys identifiers with only one redox state ----


# Functions #

url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))

# Define the theme #

theme_custom <- function(base_size = 8){
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      axis.ticks =  element_line(colour = "black"),
      panel.background = element_blank(), 
      panel.border = element_rect(fill = NA, color = "black", size = 0.5), 
      panel.grid.major = element_blank(), 
      # panel.grid.minor = element_blank(),
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      # axis.line.x = element_line(color="black", size = 0.5),
      # axis.line.y = element_line(color="black", size = 0.5)
    )
}

# Principal component analysis (PCA) #

data3 %>% plot_pca(., group) +
  theme_custom()

# Adjust dataframe for faceting #

data4 <- data3 %>% 
  left_join(., pept9, by = "Identifier1") %>% 
  select(1, matches("_P|_FDR|_FC"), HorL, RedoxState) %>% 
  filter(., RedoxState == 1) %>% 
  mutate(`REDorOX` = if_else(`HorL` == 1,
                             "Reduced",
                             "Oxidized"))

# Redefine volcano plot function #

plot_volcano <- function(data4, group, group.compare, fdr = TRUE, threshold = 2, xlimit = 10, ylimit = 8){
  # Data preparation
  temp.data <- 	data4 %>%
    #select(-unlist(group)) %>% View
    select(1, REDorOX, matches("_P|_FDR|_FC")) %>%
    gather(compare, value, -c(1,2)) %>% 
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right") %>% 
    spread(variable, value)
  
  # Check if FDR-adjustment was applied
  if (fdr == TRUE){
    temp.data <- temp.data %>%
      mutate(significance = FDR)
    
  } else {
    temp.data <- temp.data %>%
      mutate(significance = P)
    
  }
  
  # Set significance types
  temp.data <- temp.data %>%
    mutate(down = if_else(FC <= -log2(threshold) & significance < 0.05, 1, 0),
           up = if_else(FC >= log2(threshold) & significance < 0.05, 1, 0),
           type = if_else(down == 1, "down", if_else(up == 1, "up", "same")))
  
  # Build facet titles
  temp.data <- temp.data %>%
    group_by(compare) %>%
    mutate(down = sum(down), up = sum(up)) %>%
    mutate(compare_count = paste(compare, "\nDown ", sep = "", down, " / Up ", up))
  
  # Set facet order
  temp.data <- temp.data %>%
    ungroup() %>%
    mutate(compare = factor(compare, levels = names(group.compare)),
           type = factor(type, levels = c("same", "down", "up")))
  
  #
  temp.label <- temp.data %>%
    group_by(compare, compare_count) %>%
    dplyr::count() %>%
    data.frame()
  
  # Plot
  temp.data %>%
    ggplot(., aes(x = FC, y = -log10(significance), color = type)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_color_manual(values = c("same" = "grey70", "down" = "#35b779", "up" = "#31688e")) +
    coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
    xlab(expression("log"[2]*"(fold change)")) +
    ylab(if_else(fdr == TRUE,
                 #expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
                 expression("-log"[10]*"("*italic(q)*"-value)"),
                 expression("-log"[10]*"("*italic(p)*"-value)"))) +
    facet_grid(REDorOX ~ compare_count) +
    #geom_text(data = temp.label, aes(x = 0, y = Inf, label = compare_count), inherit.aes = FALSE) +
    guides(color = FALSE) +
    
    #scale_x_continuous(breaks = -xlimit:xlimit) +
    #scale_y_continuous(breaks = -ylimit:ylimit) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, color = "black") +
    geom_vline(xintercept = c(-log2(threshold), log2(threshold)), linetype = 2, size = 0.5, color = "black")
  
}


# Volcano Plot #

data4 %>%
  plot_volcano(.,
               group,
               group.compare,
               fdr = TRUE,
               threshold = 2,
               xlimit = 8,
               ylimit = 4) +
  theme_custom()

# Save volcano plot #

ggsave(filename = "TMT_HLNEM_REDandOX.pdf",
       plot = last_plot(),
       dpi = 720,
       units = "in",
       width = 3, 
       height = 3)


# -----------------------------------------------------------------------------#

# Analysis for Cys with more than redox state in vivo #

# -----------------------------------------------------------------------------#

# Oxidation percentage determination # 

# Filter for the Cys identifiers with more than redox state #
temp <- data %>% 
  left_join(., pept9, by = "Identifier1") %>% 
  filter(., RedoxState == 2)

# Sum the total intensity of a Cys, both reduced and oxidized #

temp2 <- temp %>% 
  group_by(Identifier) %>% 
  summarize_each(list(sum),2:7)
 
# Adjust column names #

colnames(temp2) <- paste("Sum", colnames(temp2))
  
temp2 <- temp2 %>% 
  select(., "Identifier" = "Sum Identifier", 2:7)
   
# Add summed intensity back to the working dataframe #

temp3 <- temp %>% 
  left_join(., temp2, by = "Identifier") %>% 
  filter(., HorL == 2) %>%
  select(1, 2:7, matches("Sum"))
  
#Define the function to determine the percent oxidized of a Cys #

per_ox <- function(df, col){
    
    y <- df
    a <- colnames(df)
    for(i in col){
      
      string <- str_extract(string = a, pattern = regex(paste0("Sum ", a[i])))
      
      string <- string[!is.na(string)]  
      
      x <- df %>% 
        transmute(!!paste0(a[i], sep = " ", "Percent Oxidized") := df[,i] / df[,paste0(string)])
      
      y <- cbind.data.frame(y, x)
      
    }
    return(y)
}

# Percent oxidized of a Cys identifier #
  
temp4 <- per_ox(temp3, col = 2:7) %>% 
  select(., 1, matches("Percent Oxidized"))

temp4[is.na(temp4)] <- 0
  
  
# Define column indeces for the percent oxidized replicates for each condition #

z <- 2:4; y <- 5:7

ox_group <- list("WT" = z , "H2O2" = y)

ox_group.compare <- list("WT-H2O2" = list(z, y))

# Rename the percent oxidized columns in a simplified "Condition-Replicate" format #

temp4 <- temp4 %>%
  rename_columns(., ox_group)

# Data cleaning and imputation #

temp5 <-temp4 %>%
  clean_min(., ox_group, nonzero = 2) %>% # ??? - nonzero should be a value at least >= half the replicates
  impute_imp4p(., ox_group)

# Hypothesis testing #

temp6 <- temp5 %>%
  calculate_ttest(., ox_group.compare) %>% # Pairwise t-test
  calculate_1anova(., ox_group)

# Define the change in oxidation function #
  
calculate_ChangeOx <- function(df, ox_group.compare) {
  for (x in ox_group.compare){
    temp.data <- df
    
    temp.meandiff <- rowMeans(temp.data[x[[2]]]) - rowMeans(temp.data[x[[1]]])
    
    # Add to dataframe
    df <- cbind(df, temp.name = temp.meandiff)
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- df %>%
      select(x[[1]][1]) %>%
      names() %>%
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    temp2 <- df %>%
      select(x[[2]][1]) %>%
      names() %>%
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    temp.name <- paste(temp1, temp2, sep = "-") %>% paste(., "ChangeOx", sep = "_")
    
    # Rename appended column with dynamic variable
    names(df)[names(df) == "temp.name"] <- temp.name
    
  }
  return(df)
}

# Define change in oxidation maximum function #

add_ChangeOx_max <- function(df){
  variable <- df %>%
    select(1) %>%
    names()
  
  temp.data <- df %>%
    select(1, contains("_ChangeOx")) %>%
    gather(condition, value, -1) %>%
    group_by(!!as.name(variable)) %>%
    summarize(`ChangeOx Max` = if_else(max(value) > abs(min(value)),
                               true = max(value),
                               false = min(value))) %>%
    left_join(df, ., by = variable)
  
}

# Determine the change in oxidation for a Cys identifier across conditions #

temp6 <- temp6 %>% 
  calculate_ChangeOx(., ox_group.compare) %>% 
  add_ChangeOx_max()
  

# Plot Cys identifiers with more than one redox state--------------

# Define volcano plot function for change in oxidation #

plot_volcano_ChangeOx <- function(temp6, ox_group, ox_group.compare, fdr = TRUE, threshold = 0.10, xlimit = 1, ylimit = 8){
  # Data preparation
  temp.data <- temp6 %>%
    #select(-unlist(group)) %>% View
    select(1, matches("_P|_FDR|_ChangeOx")) %>%
    gather(compare, value, -1) %>% 
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right") %>% 
    spread(variable, value)
  
  # Check if FDR-adjustment was applied
  if (fdr == TRUE){
    temp.data <- temp.data %>%
      mutate(significance = FDR)
    
  } else {
    temp.data <- temp.data %>%
      mutate(significance = P)
    
  }
  
  # Set significance types
  temp.data <- temp.data %>%
    mutate(down = if_else(ChangeOx <= -(threshold) & significance < 0.05, 1, 0),
           up = if_else(ChangeOx >= threshold & significance < 0.05, 1, 0),
           type = if_else(down == 1, "down", if_else(up == 1, "up", "same")))
  
  # Build facet titles
  temp.data <- temp.data %>%
    group_by(compare) %>%
    mutate(down = sum(down), up = sum(up)) %>%
    mutate(compare_count = paste(compare, "\nDown ", sep = "", down, " / Up ", up))
  
  # Set facet order
  temp.data <- temp.data %>%
    ungroup() %>%
    mutate(compare = factor(compare, levels = names(ox_group.compare)),
           type = factor(type, levels = c("same", "down", "up")))
  
  #
  temp.label <- temp.data %>%
    group_by(compare, compare_count) %>%
    dplyr::count() %>%
    data.frame()
  
  # Plot
  temp.data %>%
    ggplot(., aes(x = ChangeOx, y = -log10(significance), color = type)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_color_manual(values = c("same" = "grey70", "down" = "#35b779", "up" = "#31688e")) +
    coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
    xlab(expression("Change in Oxidation")) +
    ylab(if_else(fdr == TRUE,
                 #expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
                 expression("-log"[10]*"("*italic(q)*"-value)"),
                 expression("-log"[10]*"("*italic(p)*"-value)"))) +
    facet_grid(~ compare_count) +
    #geom_text(data = temp.label, aes(x = 0, y = Inf, label = compare_count), inherit.aes = FALSE) +
    guides(color = FALSE) +
    
    #scale_x_continuous(breaks = -xlimit:xlimit) +
    #scale_y_continuous(breaks = -ylimit:ylimit) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, color = "black") +
    geom_vline(xintercept = c(-(threshold), (threshold)), linetype = 2, size = 0.5, color = "black")
  
}

# Volcano plot for change in oxidation #

temp6 %>%
  plot_volcano_ChangeOx(.,
               ox_group,
               ox_group.compare,
               fdr = TRUE,
               threshold = .10,
               xlimit = 1,
               ylimit = 4) +
  theme_custom()
  
# Save volcano plot for change in oxidation #

ggsave(filename = "TMT_HLNEM_ChangeInOxidation.pdf",
       plot = last_plot(),
       dpi = 720,
       units = "in",
       width = 3, 
       height = 3)

# Jitter plot for change in oxidation #



plot_boxplot_ox <- function(df, group){
  temp.df <- df %>%
    select(1, group %>% flatten_int()) %>%
    gather(sample, abundance, -1)
  
  temp.df %>%  
    ggplot(., aes(x = sample, y = abundance, fill = sample)) +
    #geom_jitter(size = 1.0, height = 0, alpha = 0.8)+
    geom_boxplot(color = "black",
                 outlier.shape = NA,
                 size = 0.5,
                 alpha = 0.8) +
    scale_fill_viridis_d(limits = names(sample))+
    xlab(expression(Sample)) +
    ylab(expression(`Percent Oxidized`)) +
    guides(color = FALSE, fill = FALSE)
  
}

temp6 %>% 
  plot_boxplot_ox(., 
                  ox_group) +
  theme_custom()


ggsave(filename = "TMT_HLNEM_PercentOxidized.pdf",
       plot = last_plot(),
       dpi = 720,
       units = "in",
       width = 3, 
       height = 3)



plot_histogram_ox <- function(df, group){
  temp.df <- df %>%
    select(1, group %>% flatten_int()) %>%
    gather(sample, abundance, -1) %>% 
    separate(sample, into = c("Condition", "sample"), sep = "-") %>% 
    group_by(Identifier1, Condition) %>% 
    mutate(`Average abundance` = mean(abundance)) %>% 
    top_n(1, sample) %>% 
    dplyr::slice(1) %>% 
    ungroup()
  
  
  
  temp.df %>%  
    ggplot(., aes(x = abundance, color = Condition, fill = Condition)) +
    #geom_jitter(size = 1.0, height = 0, alpha = 0.8)+
    geom_histogram(binwidth = 0.1, position = position_dodge2(padding = 0.2, preserve = "single" )) +
    scale_color_viridis_d(begin = 0.1, end = 0.7) +
    scale_fill_viridis_d(begin = 0.1, end = 0.7, alpha = 0.8)+
    xlab(expression(`Percent Oxidized`)) +
    ylab(expression(`Total Count`))
    
}

temp6 %>% 
  plot_histogram_ox(., ox_group) +
  theme_custom() +
  theme(legend.position = "right")

ggsave(filename = "TMT_HLNEM_PercentOxidized_histogram_averages.png",
       plot = last_plot(),
       dpi = 720,
       units = "in",
       width = 6, 
       height = 3)





