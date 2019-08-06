library(skimr)
library(scales)

# Functions----
summarize_count <- function(df, col){
  df %>%
    count(!!as.name(col)) %>%
    filter(!is.na(!!as.name(col))) %>%
    mutate(freq = n / sum(n)) %>%
    arrange(desc(n))
  
}

plot_count <- function(df, col, threshold = 0){
  library(scales)
  
  temp.df <- df %>%
    count(!!as.name(col)) %>%
    filter(!is.na(!!as.name(col))) %>%
    mutate(freq = n / sum(n))
  
  # Threshold and factorize
  temp.df <- temp.df %>%
    filter(n > threshold) %>%
    mutate(term = str_trunc(!!as.name(col), 25)) %>%
    mutate(term = fct_reorder(term, n, .desc = TRUE))
  
  temp.df %>%
    ggplot(., aes(x = term, y = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.5, size = 8) +
    geom_text(aes(y = 0, label = scales::percent(freq)), vjust = 2, size = 8) +    
    scale_y_continuous(expand = c(0.1, 0.1)) +
    xlab(col) +
    ylab("Count")
  
}






# Proteins
temp.data <- data5

temp.data %>%
  group_by(Accession)


# Sites per identifier
temp.data <- temp.data %>%
  mutate("Sites_Identifier" = str_count(Sites, "-") + 1)

# Identifiers per accession
temp.data <- temp.data %>%
  add_count(Accession, name = "Identifiers_Accession")

# Site per accession
temp.data <- temp.data %>%
  separate_rows(Sites, sep = "-") %>%
  count(Accession, Sites) %>%
  count(Accession, name = "Sites_Accession") %>%
  left_join(temp.data, ., by = "Accession")


# Sites per protein
temp.data %>%
  group_by(Accession) %>%
  slice(1) %>%
  ungroup()

temp.data %>%
  summarize_count(., col = "Accession")




png("figure.png", width = 12, height = 9, units = "in", res = 300)

dev.off()






