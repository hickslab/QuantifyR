library(skimr)
library(scales)

# Functions----

summarize_cv <- function(df, group){
  # Define experiment
  variable <- df %>%
    select(1) %>%
    names()
  
  # Select abundance columns
  temp.df <- df %>%
    select(1, group %>% flatten_int())
  
  # Calculate within-condition CV per observation
  temp.df <- temp.df %>%
    gather(sample, abundance, -1) %>%
    separate(sample, into = c("condition", "replicate"), sep = "-") %>%
    group_by(!!as.name(variable), condition) %>%
    summarize(mean = mean(abundance),
              sd = sd(abundance),
              cv = (sd / mean) * 100) %>%
    mutate(class = if_else(str_detect(!!as.name(variable), "_UPS"),
                           "UPS1",
                           "Background"))
  
  temp.df %>%
    ggplot(., aes(x = mean, y = sd)) +
    #ggplot(., aes(x = log2(mean), cv)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "loess", size = 2) +
    #ylim(0, 100) +
    #scale_y_continuous(breaks = seq(0, 200, 50)) +
    #scale_color_distiller() +
    facet_grid(fct_relevel(condition, names(group)) ~ class) +
    #labs(x = expression("log"[10]*"(mean abundance)"), y = "CV (%)") +
    guides(color = FALSE) +
    theme_custom() %+replace%
    theme(strip.text.y = element_text(angle = 0))
  
  
  a <- temp.df %>%
    group_by(condition) %>%
    summarize(cv_mean = mean(cv),
              cv_median = median(cv),
              cv_90 = quantile(cv, 0.9))
  
  a
  
}


summarize_volcano <- function(df, group){
  temp.df <- 	df %>%
    #select(-unlist(group)) %>% View
    select(1, matches("_P|_FDR|_FC")) %>%
    gather(compare, value, -1) %>%
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right") %>%
    spread(variable, value) %>%
    mutate(class = if_else(str_detect(!!as.name(variable), "_UPS"), "UPS1", "Background"))
  
  temp.df %>%
    ggplot(., aes(x = FC, y = -log10(FDR))) +
    geom_point(aes(color = class), size = 5, alpha = 0.7) +
    scale_color_manual(values = c("grey70", "red")) +
    #xlim(-3, 3) +
    scale_x_continuous(limits = c(-3, 3), breaks = c(-3:3)) +
    #coord_cartesian(xlim = -3:3) +
    labs(x = expression("log"[2]*"(fold change)"), y = expression("-log"[10]*"(FDR)"), color = "") +
    facet_grid(fct_relevel(compare, names(group.compare)) ~ .) +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = 2, size = 1) +
    #guides(color = FALSE) +
    theme_custom() %+replace%
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "top") -> p
  
}


summarize_ttest_fdr <- function(df){
  temp.df <- df %>%
    select(1, matches("_FDR")) %>%
    gather(compare, value, -1) %>%
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right") %>%
    spread(variable, value) %>%
    mutate(class = if_else(str_detect(!!as.name(variable), "_UPS"), "UPS1", "Background"))
  
  temp.df %>%
    filter(FDR < 0.05) %>%
    count(compare, class) %>%
    spread(class, n)
  
}


summarize_fc <- function(df){
  # Define experiment
  variable <- df %>%
    select(1) %>%
    names()
  
  temp.df <- df %>%
    select(1, matches("_FC")) %>%
    gather(compare, value, -1) %>%
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right") %>%
    spread(variable, value) %>%
    mutate(class = if_else(str_detect(!!as.name(variable), "_UPS"), "UPS1", "Background"))
  
  temp.df %>%
    mutate(FC = 2^FC) %>%
    group_by(compare, class) %>%
    summarize(median = median(FC),
              mean = mean(FC),
              sd = sd(FC)) %>%
    ungroup()
  
  
  temp.df %>%
    mutate(FC = 2^FC) %>%
    ggplot(., aes(x = FC, fill = class)) +
    geom_histogram(bins = 100, alpha = 0.5) +
    scale_x_continuous(limits = c(-1, 5), breaks = c(-3:3)) +
    facet_grid(compare ~ .) +
    guides(fill = FALSE) +
    theme_custom()
  
}


summarize_z <- function(df, group){
  # Define experiment
  variable <- df %>%
    select(1) %>%
    names()
  
  # Select abundance columns
  temp.df <- df %>%
    select(1, group %>% flatten_int())
  
  # Calculate mean condition abundance
  temp.df <- temp.df %>%
    gather(replicate, abundance, -1) %>%
    separate(replicate, into = c("condition", "replicate"), sep = "-") %>%
    group_by(!!as.name(variable), condition) %>%
    summarize(mean = mean(abundance)) %>%
    ungroup() %>%
    spread(., condition, mean)
  
  # Z-score rowwise normalization
  temp.df[-1] <- temp.df %>%
    select(-1) %>%
    apply(., 1, scale) %>%
    t()
  
  temp.df <- temp.df %>%
    gather(condition, scaled, -1) %>%
    mutate(class = if_else(str_detect(!!as.name(variable), "_UPS"), "UPS1", "Background"))
  
  temp.df %>%
    ggplot(., aes(x = fct_relevel(condition, names(group)), y = scaled)) +
    geom_line(aes(group = !!as.name(variable), color = class), alpha = 0.5, size = 1) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
    scale_color_manual(values = c("grey70", "red")) +
    facet_wrap(~ class) +
    labs(x = "Condition", y = "Z-score") +
    guides(color = FALSE) +
    theme_custom() -> p
  
  temp.df %>%
    ggplot(., aes(x = fct_relevel(condition, names(group)), y = scaled)) +
    geom_jitter(aes(color = class), height = 0, size = 3, alpha = 0.1) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
    scale_color_manual(values = c("grey70", "red")) +
    guides(color = FALSE, fill = FALSE) +
    facet_wrap(~ class) +
    labs(x = "Condition", y = "Z-score") +
    theme_custom()
  
}



summarize_count <- function(df, col){
  df %>%
    count(!!as.name(col)) %>%
    filter(!is.na(!!as.name(col))) %>%
    mutate(freq = n / sum(n)) %>%
    arrange(desc(n))
  
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







