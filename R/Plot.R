theme_custom <- function(base_size = 32){
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      axis.ticks =  element_line(colour = "black"),
      panel.background = element_blank(), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}


plot_jitter <- function(df, group){
  temp.df <- df %>%
    select(1, group %>% flatten_int()) %>%
    gather(sample, abundance, -1)
  
  temp.df %>%  
    ggplot(., aes(x = sample, y = abundance, color = sample)) +
    geom_jitter(height = 0, alpha = 0.5) +
    geom_boxplot(color = "black",
                 fill = NA,
                 outlier.shape = NA,
                 size = 1.5) +
    guides(color = FALSE, fill = FALSE)
  
}


plot_corr <- function(df, group){
  temp.df <- df %>%
    select(1, group %>% flatten_int()) %>%
    select(-1) %>%
    cor() %>%
    as_tibble(rownames = "samples") %>%
    gather(samples2, r, -1)
  
  temp.df %>%
    ggplot(., aes(x = samples, y = samples2)) +
    geom_raster(aes(fill = r)) +
    scale_fill_distiller(palette = "Reds", limits = c(0.8, 1)) +
    geom_text(aes(label = round(r, 2)), size = 6) +
    guides(fill = FALSE) +
    labs(x = NULL, y = NULL)
  
}


plot_dendrogram <- function(df, group, k = 3){
  # Load packages
  library(dendextend)
  
  temp.df <- df %>%
    select(group %>% flatten_int()) %>%
    scale() %>%
    t() %>%
    dist() %>%
    hclust() %>%
    #cutree(., 3) %>%
    as.dendrogram()
  
  temp.df %>%
    set("branches_k_color", k = 3) %>%
    plot(xlab = "Distance", ylab = "Replicate")
  
  temp.df %>%
    rect.dendrogram(k)
  
  
  library(ggdendro)
  
  temp.df <- df %>%
    select(group %>% flatten_int()) %>%
    scale() %>%
    t() %>%
    dist() %>%
    hclust()
  
  temp.df %>%
    as.dendrogram() %>%
    dendro_data() %>%
    ggdendrogram(rotate = TRUE, size = 2) +
    labs(x = "", y = "Distance")
    theme_custom()
  
  
    
  
}


plot_pca <- function(df, group){
  temp.df <- df %>%
    select(group %>% flatten_int())
  
  temp.df %>%
    prcomp() %>%
    .$rotation %>%
    data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, into = c("condition", "replicate"), sep = "-") %>%
    mutate(condition = fct_relevel(condition, names(group))) %>%
    
    ggplot(., aes(x = PC1, y = PC2)) +
    geom_point(aes(color = condition, shape = condition), alpha = 0.5, size = 18) +
    geom_text(aes(label = replicate), color = "black", size = 10) +
    #stat_ellipse(aes(color = condition)) +
    #scale_color_discrete(limits = names(group)) +
    labs(color = "Condition", shape = "Condition")
  
}


plot_volcano <- function(data3, group, group.compare, fdr = TRUE, threshold = 2, xlimit = 10, ylimit = 8){
  # Data preparation
  temp.data <- 	data3 %>%
	  #select(-unlist(group)) %>% View
	  select(1, matches("_P|_FDR|_FC")) %>%
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
	  geom_point(size = 5, alpha = 0.8) +
	  scale_color_manual(values = c("same" = "grey70", "down" = "blue", "up" = "red")) +
	  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
	  xlab(expression("log"[2]*"(fold change)")) +
	  ylab(if_else(fdr == TRUE,
	               #expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
	               expression("-log"[10]*"("*italic(q)*"-value)"),
	               expression("-log"[10]*"("*italic(p)*"-value)"))) +
	  facet_wrap(~ compare_count) +
	  #geom_text(data = temp.label, aes(x = 0, y = Inf, label = compare_count), inherit.aes = FALSE) +
	  guides(color = FALSE) +
	  
	  #scale_x_continuous(breaks = -xlimit:xlimit) +
	  #scale_y_continuous(breaks = -ylimit:ylimit) +
	  
	  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 1, color = "black") +
	  geom_vline(xintercept = c(-log2(threshold), log2(threshold)), linetype = 2, size = 1, color = "black")

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


plot_box <- function(df){
  temp.df %>%
    gather(sample, abundance, -1) %>%
    separate(sample, into = c("condition", "replicate"), sep = "-", remove = FALSE) %>%
    mutate(sample = fct_relevel(sample, names(group))) %>%
    
    ggplot(., aes(x = sample, y = abundance, color = condition)) +
    geom_jitter(alpha = 0.5, height = 0) +
    geom_boxplot(color = "black",
                 fill = NA,
                 outlier.shape = NA,
                 size = 1.5) +
    guides(color = FALSE, fill = FALSE) +
    #coord_flip() +
    theme_custom() +
    labs(x = "Replicate", y = expression("log"[2]*"(abundance)"))
  
}


plot_hclust <- function(df, group, k = 3, type = "jitter"){
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
    select(-1) %>%
    
    mutate(cluster = dist(.) %>%
             hclust() %>%
             cutree(k)) %>%
    
    mutate(cluster = LETTERS[cluster]) %>%
    
    group_by(cluster) %>%
    mutate(cluster_count = paste(cluster, "\n", sep = "", n())) %>%
    ungroup() %>%
    
    rownames_to_column() %>%
    
    gather(condition, scaled, -rowname, -cluster, -cluster_count)
  
  # Set condition order
  temp.df <- temp.df %>%
    mutate(condition = factor(condition, level = names(group)))
  
  
  temp.df <- temp.df %>%
    mutate(breaks = group_indices(., condition))
    
  # Plot
  p <- temp.df %>%
    ggplot(., aes(x = factor(breaks), y = scaled)) +
    geom_jitter(aes(color = condition), height = 0) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
    geom_smooth(aes(x = jitter(breaks)), method = "loess", size = 2) +
    scale_x_discrete(breaks = 1:length(group), labels = names(group)) +
    guides(color = FALSE, fill = FALSE) +
    facet_wrap(~ cluster_count) +
    labs(x = "Condition", y = "Z-score")
  
  # Line type
  if (type == "line"){
    p <- temp.df %>%
      ggplot(., aes(x = factor(breaks), y = scaled)) +
      geom_line(aes(group = rowname, color = cluster_count), alpha = 0.5, size = 1) +
      geom_smooth(aes(x = jitter(breaks)), method = "loess", size = 1.5) +
      geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
      scale_x_discrete(breaks = 1:length(group), labels = names(group)) +
      guides(color = FALSE, fill = FALSE) +
      facet_wrap(~ cluster_count) +
      labs(x = "Condition", y = "Z-score")
    
  }
  return(p)
  
}

 
plot_heatmap <- function(df){

  df %>%
    ggplot(., aes(variable, reorder(Identifier, -value))) +
    #ggplot(., aes(variable, Identifier)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("")
  
}


plot_GO <- function(df, column = "Gene ontology", top = 5){
  # Select GO columns
  temp.df <- df %>%
    select(Accession, contains(column)) %>%
    distinct(., Accession, .keep_all = TRUE) %>%
    gather(column, term, -1) %>%
    
    mutate(column = str_split(column, "\\(", simplify = TRUE)[, 2] %>% str_sub(., end = -2)) %>%
    
    filter(!is.na(term)) %>%
    separate_rows(term, sep = "; ") %>%
    mutate(term = str_sub(term, end = -14)) %>%
    mutate(term = str_trunc(term, 50))
    
  # Count terms
  temp.df <- temp.df %>%
    count(column, term)
  
  # Filter for top terms
  temp.df <- temp.df %>%
    group_by(column) %>%
    top_n(., top, n) %>%
    arrange(desc(n)) %>%
    dplyr::slice(1:top)
  
  # Plot
  temp.df %>%
    mutate(Cluster = "Proteins") %>%
    ggplot(., aes(x = reorder(term, n), y = Cluster)) +
    geom_raster(aes(fill = column, alpha = n)) +
    geom_text(aes(label = n), size = 10) +
    scale_alpha_continuous(range = c(0.5, 1)) +
    facet_grid(column ~ ., scales = "free") +
    guides(color = FALSE, fill = FALSE, alpha = FALSE) +
    labs(x = NULL, y = NULL) +
    coord_flip() +
    theme_custom()
  
}


plot_GO_cluster <- function(df, column = "Gene ontology", top = 5){
  # Select columns
  temp.df <- df %>%
    select(Accession, Cluster, contains(column)) %>%
    distinct(., Accession, .keep_all = TRUE) %>%
    gather(column, term, -1, -Cluster) %>%
    
    mutate(column = str_split(column, "\\(", simplify = TRUE)[, 2] %>% str_sub(., end = -2)) %>%
    
    filter(!is.na(term)) %>%
    separate_rows(term, sep = "; ") %>%
    mutate(term = str_sub(term, end = -14)) %>%
    mutate(term = str_trunc(term, 50))
  
  # Count terms
  temp.df <- temp.df %>%
    count(column, term, Cluster)
  
  # Filter for top terms in each cluster
  temp.df <- temp.df %>%
    group_by(column, Cluster) %>%
    top_n(., top, n) %>%
    dplyr::slice(1:top) %>%
    mutate(rank = TRUE) %>%
    ungroup() %>%
    select(term, rank) %>%
    left_join(temp.df, ., by = "term") %>%
    filter(rank == TRUE)
  
  # Plot
  temp.df %>%
    ggplot(., aes(x = reorder(term, n), y = Cluster)) +
    geom_raster(aes(fill = column, alpha = n)) + geom_text(aes(label = n), size = 6) + scale_alpha_continuous(range = c(0.5, 1)) +
    #geom_point(size = 10, color = "grey90") + geom_point(aes(color = Cluster, size = n)) + scale_size_continuous(range = c(2, 10)) +
    facet_grid(column ~ ., scales = "free") +
    guides(color = FALSE, fill = FALSE, alpha = FALSE, size = FALSE) +
    labs(x = NULL, y = "Cluster", alpha = "Proteins") +
    coord_flip()
  
}


plot_GO_hclust <- function(df, group, column = "Gene ontology (biological process)", threshold = 10){
  # Define experiment
  variable <- df %>%
    select(1) %>%
    names()
  
  # Select abundance columns
  temp.data <- df %>%
    select(1, group %>% flatten_int())
  
  # Calculate mean condition abundance
  temp.data <- temp.data %>%
    gather(replicate, abundance, -1) %>%
    separate(replicate, into = c("condition", "replicate"), sep = "-") %>%
    group_by(!!as.name(variable), condition) %>%
    summarize(mean = mean(abundance)) %>%
    ungroup() %>%
    spread(., condition, mean)
  
  # Z-score rowwise normalization
  temp.data[-1] <- temp.data %>%
    select(-1) %>%
    apply(., 1, scale) %>%
    t()

  # Add annotation column
  temp.data <- df %>%
    select(1, column) %>%
    dplyr::rename(column = !!as.name(column)) %>%
    left_join(temp.data, ., by = variable)
  
  # Melt by annotation
  temp.data <- temp.data %>%
    mutate(column = str_sub(column, end = -2)) %>%
    separate_rows(column, sep = ";") %>%
    gather(condition, scaled, -1, -column) %>%
    filter(!is.na(column))
  
  # Count and threshold
  temp.data <- temp.data %>%
    count(!!as.name(variable), column) %>%
    count(column) %>%
    left_join(temp.data, ., by = "column") %>%
    mutate(column_n = str_c(column, "\n", n)) %>%
    filter(n > threshold)
    
  # Plot
  temp.data %>%
    ggplot(., aes(x = factor(condition), y = scaled)) +
    geom_line(aes(group = !!as.name(variable), color = column), size = 1.5) +
    #geom_smooth(aes(x = jitter(breaks)), method = "loess", size = 1.5) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
    #scale_x_discrete(breaks = 1:length(group), labels = names(group)) +
    guides(color = FALSE, fill = FALSE) +
    facet_wrap(~ fct_reorder(column_n, n, .desc = TRUE)) +
    labs(x = "Condition", y = "Z-score")
  
}


plot_GO_FC <- function(df, group, column = "Gene ontology (biological process)", threshold = 5){
  # Define experiment
  variable <- df %>%
    select(1) %>%
    names()
  
  # Select data
  temp.df <- df %>%
    select(1, contains("_FC"))
           
  temp.df <- df %>%
    select(1, column) %>%
    dplyr::rename(column = !!as.name(column)) %>%
    left_join(temp.df, ., by = variable)
  
  # Melt by annotation
  temp.df <- temp.df %>%
    #mutate(column = str_sub(column, end = -2)) %>%
    separate_rows(column, sep = "; ") %>%
    gather(condition, scaled, -1, -column) %>%
    filter(!is.na(column))
  
  # Count and threshold
  temp.df <- temp.df %>%
    count(!!as.name(variable), column) %>%
    count(column) %>%
    left_join(temp.df, ., by = "column") %>%
    mutate(column_n = str_c(column, "\n", n)) %>%
    filter(n > threshold)
  
  # Plot
  temp.df %>%
    ggplot(., aes(x = factor(condition), y = scaled)) +
    geom_line(aes(group = !!as.name(variable), color = column), size = 1.5) +
    #geom_smooth(aes(x = jitter(breaks)), method = "loess", size = 1.5) +
    geom_boxplot(color = "black", fill = NA, outlier.shape = NA, size = 1.5) +
    #scale_x_discrete(breaks = 1:length(group), labels = names(group)) +
    guides(color = FALSE, fill = FALSE) +
    facet_wrap(~ fct_reorder(column_n, n, .desc = TRUE)) +
    labs(x = "Condition", y = "Fold change")
  
}


plot_GO_heatmap <- function(., group, column = "Gene ontology", threshold = 5){
  # Select data
  temp.df <- df %>%
    select(1, contains("scaled"))
  
  temp.df <- df %>%
    select(1, column) %>%
    dplyr::rename(column = !!as.name(column)) %>%
    left_join(temp.df, ., by = variable)
  
  # Melt by annotation
  temp.df <- temp.df %>%
    #mutate(column = str_sub(column, end = -2)) %>%
    separate_rows(column, sep = "; ") %>%
    gather(condition, scaled, -1, -column) %>%
    filter(!is.na(column))
  
  # Count and threshold
  temp.df <- temp.df %>%
    count(!!as.name(variable), column) %>%
    count(column) %>%
    left_join(temp.df, ., by = "column") %>%
    mutate(column_n = str_c(column, "\n", n)) %>%
    filter(n > threshold)
  
  # Plot
  temp.df %>%
    ggplot(., aes(x = condition, y = Identifier, size = scaled)) +
    #geom_tile(aes(fill = scaled)) +
    geom_count() +
    facet_grid(fct_reorder(column_n, n, .desc = TRUE) ~ ., scales = "free") +
    guides(color = FALSE)# + theme_custom()
  
}


plot_save <- function(p, filename = "figure.png", w = 12, h = 10, dpi = 300){
  # Single column: width = 12, height = 10
  # Double colum:  width = 10, height = 5 
  
  # Initialize image file
  png(filename, width = w, height = h, units = "in", res = dpi)
  
  # Write plot to file
  print(p)
  
  # Close file
  dev.off()
  
}
