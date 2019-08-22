simulate_data <- function(data = 1000, p = 100, fold = c(2, 4), num.rep = 3, cv = 10){
  # Initiate RNG state
  set.seed(123)
  
  temp.df <- tibble(Accession = c(str_c("N", c(1:(data - p)), sep = "_"),
                                  str_c("P", c(1:p), sep = "_")))
  
  temp.df <- temp.df %>%
    mutate(Range = rnorm(n = n(), mean = 9, sd = 2))
  
  fold2 <- c(1, fold)
  
  for (x in 1:length(fold2)){
    for (y in 1:num.rep){
      temp.name <- LETTERS[x] %>%
        paste(., y, sep = "-")
      
      temp.df <- temp.df %>%
        mutate(!!temp.name := if_else(str_detect(Accession, "P_"),
                                      rnorm(n = Range, mean = Range + log2(fold2[x]), sd = (1.2^(-Range)) * (cv / 5)),
                                      rnorm(n = Range, mean = Range, sd = (1.2^(-Range) * (cv / 5)))))
    }
  }
  temp.df <- temp.df %>%
    select(-Range)

}


simulate_group <- function(num.cond = 3, num.rep = 3){
  group <- list()
  
  idx <- 2
  
  for (x in 1:num.cond){
    assign(LETTERS[x], c(idx:(idx + num.rep - 1)))
    
    group[[LETTERS[x]]] <- get(LETTERS[x])
    
    idx <- idx + num.rep
  }
  return(group)
  
}


fold = seq(2, 220, 20)
2^rnorm(11, sd = 10)
c(2, 1/2, 2, 1/2, 1)


data2 <- simulate_data(data = 1000, p = 50, num.rep = 4, fold = c(2, 4), cv = 10)
group <- simulate_group(num.cond = 3, num.rep = 4)

data2 %>%
  ggplot(., aes(x = `A-1`, y = `A-2`)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(x = c(0, 20), y = c(0, 20))


data2 %>%
  gather(condition, abundance, -1) %>%
  separate(condition, into = c("condition", "replicate"), sep = "-") %>%
  mutate(abundance = 2^abundance) %>%
  
  group_by(Accession, condition) %>%
  summarize(cv = (sd(abundance) / mean(abundance)) * 100) %>%
  
  mutate(class = if_else(str_detect(Accession, "P_"), "P", "N")) %>%
  group_by(condition, class) %>%
  summarize(median = median(cv))

  
  
temp.data2 %>%
  ggplot(., aes(`20-1`, y = `20-2`)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(x = c(0, 20), y = c(0, 20))
  

temp.data2 %>%
  gather(condition, abundance, -1) %>%
  separate(condition, into = c("condition", "replicate"), sep = "-") %>%
  mutate(abundance = 2^abundance) %>%
  
  group_by(Identifier, condition) %>%
  summarize(cv = (sd(abundance) / mean(abundance)) * 100) %>%
  group_by(condition) %>%
  summarize(median = median(cv))



# Analyze ----


# Packages
library(imp4p)
library(broom)


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Analyze.R"))


# Workflow
# Define the column indeces for replicates in each condition.
a <- 2:5; b <- 6:9; c <- 10:13
#group <- list("A" = a, "B" = b) # ???
group.compare <- list("A-B" = list(a, b),
                      "A-C" = list(a, c),
                      "B-C" = list(b, c)) # ???


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
  #filter(FDR < 0.05) %>%
  calculate_hclust(., group, k = 2) %>%
  left_join(data3, ., by = names(data3)[1])


# Plot ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/QuantifyR/master/"
source_url(paste0(url, "R/Plot.R"))


# PCA
data3 %>%
  plot_pca(., group) +
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
#data3[group %>% flatten_int()] <- 2^data3[group %>% flatten_int()]

data3 %>%
  filter(P < 0.05) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()
