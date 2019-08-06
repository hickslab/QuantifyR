rename_columns <- function(df, group){
  variable <- df %>%
    select(1) %>%
    names()
  
  names(df) <- group %>%
    names() %>%
    rep(., group %>% map(., length)) %>%
    paste(., group %>% map(., seq) %>% flatten_int(), sep = "-") %>%
    c(variable, .)
  
  return(df)
  
}


clean_min <- function(data, group = group, nonzero = 3){
  # Create vector to store index
  idx <- c()
  i <- 1
  
  # Iterate by row
  for (x in 1:nrow(data)){
    row <- data[x,]
    
    # Initiate condition
    keep <- FALSE
    
    # Iterate by replicates for each sample
    for (y in group){
      row2 <- row[, y]
      
      # Threshold how many columns can be not equal to 0
      if (sum(row2 != 0) >= nonzero){
        keep <- TRUE
        
      }
    }
    # Check on condition state
    if(keep == TRUE){
      idx[i] <- rownames(row)
      
    }
    i <- i + 1
    
  }
  return(data[which(rownames(data) %in% idx), ])
  
}


transform_data <- function(data, group = group, method = "log2"){
  if (method == "log2"){
    data[, unlist(group)] <- log2(data[, unlist(group)])
    data[, unlist(group)][data[, unlist(group)] == -Inf] <- 0
    
  } else if (method == "asinh"){
    data[, unlist(group)] <- asinh(data[, unlist(group)])
    
  }
  return(data)
  
}


normalize_data <- function(data){
  # Would this help following data reduction/filtering?
  
}


impute_imp4p <- function(data, group = group){
  # Set RNG state
  set.seed(123)
  
  # Define condition and sample number based on group
  num.cond <- length(group)
  num.sample <- length(unlist(group)) / length(group)
  
  # Create vector showing replicates for each sample
  membership <- gen.cond(num.cond, num.sample)
  
  # Store raw data to manipulate
  temp.data <- data.matrix(data[, unlist(group)])
  
  # All 0 to NA for 'imp4p' requirement
  temp.data[temp.data == 0] <- NA
  
  # Impute rows with with nonzeros in every condition
  temp.data <- impute.rand(temp.data, membership)
  
  # Impute rows with only zeros in at least one condition
  temp.data <- impute.pa(temp.data, membership, q.norm = 0)

  # Write onto original data
  data[, unlist(group)] <- data.frame(temp.data$tab.imp)
  
  return(data)
  
}


impute_norm <- function(df, group){
  # Check data type
  variable <- df %>%
    select(1) %>%
    names()
  
  
  df %>%
    gather(condition, abundance, -variable)
  
}


calculate_ttest <- function(data2, group.compare, fdr = TRUE){
  # @TODO
  
  # Dataframe with only variable column to build upon
  #temp.data <- data %>%
  #  select(1)
  
  #for (x in group.compare){
  #data %>%
  #  select(unlist(x)) %>%
  #  rownames_to_column() %>%
  #  gather(sample, abundance, -rowname) %>%
  #  mutate(cond = str_sub(sample, start = 1, end = 1)) %>%
  #  group_by(rowname) %>%
  #  mutate(p = t.test(abundance, alternative = "two.sided", var.equal = TRUE)$p.value)
  
  # @TODO: allow FDR toggle on/off
  
  
  for (x in group.compare){
    # Row-wise t-test
    p.ttest <- c()
    for (i in 1:nrow(data2)){
      row <- data2[i, ]
      
      ttest.out <- t.test(row[, x[[1]]],
                          row[, x[[2]]],
                          alternative = "two.sided",
                          var.equal = TRUE)
      
      p.ttest[i] <- ttest.out$p.value
      
    }
    
    if (fdr == TRUE){
      # Benjamini-Hochberg FDR correction
      p.ttest <- p.adjust(p.ttest, method = "BH", n = length(p.ttest))
      
    }
    # Add to dataframe
    data2 <- cbind(data2, temp.name = p.ttest)
    
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- data2 %>%
      select(x[[1]][1]) %>%
      names() %>%
      #str_sub(., start = 1, end = 1)
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
      
    temp2 <- data2 %>%
      select(x[[2]][1]) %>%
      names() %>%
      #str_sub(., start = 1, end = 1)
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    if (fdr == TRUE){
      #temp.name <- paste(temp1, temp2, "_FDR", sep = "-")
      temp.name <- paste(temp1, temp2, sep = "-") %>% paste(., "FDR", sep = "_")
      
    } else {
      #temp.name <- paste(temp1, temp2, "_P", sep = "-")
      temp.name <- paste(temp1, temp2, sep = "-") %>% paste(., "P", sep = "_")
      
    }
    
    # Rename appended column with dynamic variable
    names(data2)[names(data2) == "temp.name"] <- temp.name
    
  }
  return(data2)
  
}


calculate_ttest2 <- function(df, group.compare){
  # Check data type
  variable <- df %>%
    select(1) %>%
    names()
  
  # Initiate output dataframe
  temp.ttest <- df
  
  # Loop over pairwise comparisons
  for (i in 1:length(group.compare)){
    # Format
    temp.data <- df %>%
      select(1, group.compare[[i]] %>% flatten_int()) %>%
      gather(condition, abundance, -1) %>%
      separate(condition, into = c("condition", "replicate"), sep = "-") %>%
      group_by_(variable)
    
    # Nest and test
    temp.data <- temp.data %>%
      nest() %>%
      mutate(ttest = map(data, ~ t.test(abundance ~ condition, data = .x, 
                                        alternative = "two.sided",
                                        var.equal = TRUE)),
             summary = map(ttest, tidy))
    
    # Unnest and FDR adjust
    temp.data <- temp.data %>%
      unnest(summary) %>%
      mutate(fdr = p.adjust(p.value, method = "BH", n = length(p.value)))
    
    # Select model output
    temp.data <- temp.data %>%
      select(1, p.value, fdr)
    
    # Rename columns
    temp.name <- group.compare[i] %>% names()
    temp.data <- temp.data %>% rename_at("p.value", ~ paste(temp.name, "_P", sep = ""))
    temp.data <- temp.data %>% rename_at("fdr", ~ paste(temp.name, "_FDR", sep = ""))
    
    # Join to data
    temp.ttest <- temp.data %>%
      left_join(temp.ttest, ., by = variable)
    
  }
  # Exit
  return(temp.ttest)
  
}


calculate_1anova <- function(data2){
  # Check data type
  variable <- data2 %>%
    select(1) %>%
    names()
  
  # Format
  temp.data <- data2 %>%
    gather(condition, abundance, -1) %>%
    separate(condition, into = c("condition", "replicate"), sep = "-") %>%
    group_by(!!as.name(variable))
  
  # Nest and test
  temp.data <- temp.data %>%
    nest() %>%
    mutate(aov = map(data, ~ aov(abundance ~ condition, data = .x)),
           summary = map(aov, tidy))
  
  # Unnest and FDR adjust
  temp.data <- temp.data %>%
    unnest(summary) %>%
    filter(term == "condition") %>%
    dplyr::rename(., P = p.value) %>%
    mutate(FDR = p.adjust(P, method = "BH", n = length(P)))
  
  # Join to data
  temp.data <- temp.data %>%
    select(1, P, FDR) %>%
    left_join(data2, ., by = variable)
  
  # Exit
  return(temp.data)
  
}

calculate_fc <- function(data, group.compare, difference = TRUE){
  for (x in group.compare){
    temp.data <- data
    
    # Specify calculation based on 'difference' argument
    if (difference == TRUE){
      temp.fc <- rowMeans(temp.data[x[[2]]]) - rowMeans(temp.data[x[[1]]])
      
    } else if (difference == FALSE){
      temp.fc <- 2^rowMeans(temp.data[x[[2]]]) / 2^rowMeans(temp.data[x[[1]]])
      
    }
    # Add to dataframe
    data <- cbind(data, temp.name = temp.fc)
    
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- data %>%
      select(x[[1]][1]) %>%
      names() %>%
      #str_sub(., start = 1, end = 1)
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    temp2 <- data %>%
      select(x[[2]][1]) %>%
      names() %>%
      #str_sub(., start = 1, end = 1)
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    #temp.name <- paste(temp1, temp2, "_FC", sep = "")
    temp.name <- paste(temp1, temp2, sep = "-") %>% paste(., "FC", sep = "_")
    
    
    # Rename appended column with dynamic variable
    names(data)[names(data) == "temp.name"] <- temp.name
    
  }
  return(data)
  
}


add_fc_max <- function(df){
  variable <- df %>%
    select(1) %>%
    names()
  
  temp.data <- df %>%
    select(1, contains("_FC")) %>%
    gather(condition, value, -1) %>%
    group_by(!!as.name(variable)) %>%
    summarize(FC_max = if_else(max(value) > abs(min(value)),
                                   true = max(value),
                                   false = min(value))) %>%
    left_join(df, ., by = variable)
  
}



calculate_hclust <- function(df, group, k = 3){
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
  
  # Rename Z-score columns
  names(temp.data)[-1] <- names(temp.data)[-1] %>% paste0(., "_Z")
  
  temp.data <- temp.data %>%
    select(-1) %>%
    
    mutate(Cluster = dist(.) %>%
             hclust() %>%
             cutree(k)) %>%
    
    mutate(Cluster = LETTERS[Cluster]) %>%
    
    bind_cols(temp.data[1], .)
  
  # Exit
  return(temp.data)
  
}
