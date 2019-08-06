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


plot_corr <- function(data, group = group) {
  # Load packages
  library(broom)
  
  data2 %>%
    select(-1) %>%
    cor() %>%
    tidy() %>%
    gather(replicate, cor, -.rownames) %>%
    ggplot(., aes(x = .rownames, y = replicate, fill = cor)) +
    geom_tile()
  
}


plot_dendrogram <- function(df, group, k = 3){
  # Load packages
  library(dendextend)
  
  temp.data <- df %>%
    select(group %>% flatten_int()) %>%
    scale() %>%
    t() %>%
    dist() %>%
    hclust() %>%
    #cutree(., 3) %>%
    as.dendrogram()
  
  temp.data %>%
    set("branches_k_color", k = 3) %>%
    plot(xlab = "Distance", ylab = "Replicate")
  
  temp.data %>%
    rect.dendrogram(k)
  
}


plot_pca <- function(df, group){
  temp.data <- df %>%
    select(group %>% flatten_int())
  
  temp.data %>%
    prcomp() %>%
    .$rotation %>%
    data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, into = c("condition", "replicate"), sep = "-") %>%
    
    ggplot(., aes(x = PC1, y = PC2, color = condition, label = replicate)) +
    geom_point(alpha = 0.5, size = 20) +
    geom_text(color = "black", size = 10) +
    scale_color_discrete(limits = names(group)) +
    labs(color = "Condition")
  
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
    count() %>%
    data.frame()
  
	# Plot
	temp.data %>%
	  ggplot(., aes(x = FC, y = -log10(significance), color = type)) +
	  geom_point(size = 3, alpha = 0.8) +
	  scale_color_manual(values = c("same" = "grey70", "down" = "blue", "up" = "red")) +
	  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
	  xlab(expression("log"[2]*"(fold change)")) +
	  ylab(if_else(fdr == TRUE,
	               #expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
	               expression( "-log"[10]*"(FDR)"),
	               expression( "-log"[10]*"("*italic(p)*"-value)"))) +
	  facet_wrap(~ compare_count) +
	  #geom_text(data = temp.label, aes(x = 0, y = Inf, label = compare_count), inherit.aes = FALSE) +
	  guides(color = FALSE) +
	  
	  #scale_x_continuous(breaks = -xlimit:xlimit) +
	  scale_y_continuous(breaks = -ylimit:ylimit) +
	  
	  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 1, color = "black") +
	  geom_vline(xintercept = c(-log2(threshold), log2(threshold)), linetype = 2, size = 1, color = "black")

}


plot_hclust <- function(df, group, k = 3, type = "jitter"){
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
  
  temp.data <- temp.data %>%
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
  temp.data <- temp.data %>%
    mutate(condition = factor(condition, level = names(group)))
  
  
  temp.data <- temp.data %>%
    mutate(breaks = group_indices(., condition))
    #mutate(breaks = condition %>% as.numeric())
    
  # Plot
  p <- temp.data %>%
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
    p <- temp.data %>%
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


plot_GO <- function(df, top = 5){
  # Select GO columns
  temp.df <- df %>%
    #select(contains("Gene ontology")) %>%
    #gather(column, term) %>%
    
    select(Accession, contains("Gene ontology")) %>%
    distinct(., Accession, .keep_all = TRUE) %>%
    gather(column, term, -1) %>%
    
    mutate(column = str_split(column, "\\(", simplify = TRUE)[, 2] %>% str_sub(., end = -2)) %>%
    
    filter(!is.na(term)) %>%
    separate_rows(term, sep = "; ") %>%
    mutate(term = str_sub(term, end = -14)) %>%
    mutate(term = str_trunc(term, 50))
    
  # Count terms
  temp.df %>%
    count(column, term) %>%

    group_by(column) %>%
    top_n(., top, n) %>%
    dplyr::slice(1:top) %>%

    ggplot(., aes(x = reorder(term, n), y = n, fill = column)) +
    geom_bar(stat = "identity") +
    facet_grid(column ~ ., scales = "free") +
    guides(fill = FALSE) +
    labs(x = NULL, y = "Proteins") +
    coord_flip()
  
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


plot_box <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(variable, value, fill = group)) +
    geom_boxplot(outlier.alpha = 0.1) +
    facet_wrap(~group, scales = "free_x") +
    scale_fill_discrete()
  
}


plot_violin <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(variable, value, fill = group)) +
    geom_violin(trim = TRUE) +
    #geom_jitter(alpha = 0.1) +
    scale_fill_discrete()
  
}


plot_joy <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(value, variable, height = ..density.., fill = group)) +
    geom_joy(stat = "density", scale = 3) +
    #theme_joy(base_size = 14) +
    theme(axis.text = element_text(vjust = 0, color = "black")) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    xlab(expression("log"[2]*"(Abundance)")) +
    #facet_wrap(~group, scales = "free_y") +
    ylab("")
  
}


plot_fc_density <- function(df, group){
  #df: calculate_fc()
  
  # Data preparation
  temp.data <- 	df %>%
    select(1, contains("_FC")) %>%
    gather(compare, value, -1) %>%
    separate(compare,
             sep = "_",
             into = c("compare", "variable"),
             extra = "merge",
             fill = "right")
  
  ggplot(., aes(x = value, fill = compare)) +
    geom_density(alpha = 0.5)
  
  
}


plot_2venn <- function(df){
  library(ggforce)
  
  df <- tibble(label = c("A", "B"),
               x = c(-1, 1),
               y = c(0, 0))
  
  df %>%
    ggplot(., aes(x0 = x, y0 = y, r = 2, fill = label)) +
    geom_circle(alpha = 0.5)
  
}


plot_save <- function(p, dpi = 300){
  png("figure.png", width = 12, height = 9, units = "in", res = dpi)
  
  print(p)
  
  dev.off()
  
}
