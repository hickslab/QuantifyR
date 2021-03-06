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


simulate_compare <- function(num.cond = 3, num.rep = 3){
  group.compare <- list()
  
  assign(LETTERS[1], c(2:(2 + num.rep - 1)))
  
  idx <- 2 + num.rep
  
  for (x in 1:num.cond){
    assign(LETTERS[x + 1], c(idx:(idx + num.rep - 1)))
    
    group.compare[[paste(LETTERS[1], LETTERS[x + 1], sep = "-")]] <- list(get(LETTERS[1]), get(LETTERS[x + 1]))
    
    idx <- idx + num.rep
  }
  return(group.compare)
  
}


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


transform_data <- function(data, group, method = "log2"){
  if (method == "log2"){
    data[, unlist(group)] <- log2(data[, unlist(group)])
    data[, unlist(group)][data[, unlist(group)] == -Inf] <- 0
    
  } else if (method == "asinh"){
    data[, unlist(group)] <- asinh(data[, unlist(group)])
    
  }
  return(data)
  
}


impute_imp4p <- function(data, group){
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


calculate_ttest <- function(df, group.compare){
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
      group_by(!!as.name(variable))
    
    # Nest and test
    temp.data <- temp.data %>%
      nest() %>%
      mutate(ttest = map(data, ~ t.test(abundance ~ condition, data = .x, 
                                        alternative = "two.sided",
                                        var.equal = TRUE)),
             summary = map(ttest, tidy)) %>%
      unnest(summary)
    
    # FDR adjust
    temp.data2 <- temp.data %>%
      ungroup() %>%
      mutate(fdr = p.adjust(p.value, method = "BH", n = length(p.value)))
    
    # Select model output
    temp.data3 <- temp.data2 %>%
      select(1, p.value, fdr)
    
    # Rename columns
    temp.name <- group.compare[i] %>% names()
    temp.data3 <- temp.data3 %>% rename_at("p.value", ~ paste(temp.name, "_P", sep = ""))
    temp.data3 <- temp.data3 %>% rename_at("fdr", ~ paste(temp.name, "_FDR", sep = ""))
    
    # Join to data
    temp.ttest <- temp.data3 %>%
      left_join(temp.ttest, ., by = variable)
    
  }
  # Exit
  return(temp.ttest)
  
}


calculate_1anova <- function(df, group){
  # Check data type
  variable <- df %>%
    select(1) %>%
    names()
  
  # Format
  temp.df <- df %>%
    select(1, group %>% flatten_int()) %>%
    gather(condition, abundance, -1) %>%
    separate(condition, into = c("condition", "replicate"), sep = "-") %>%
    group_by(!!as.name(variable))
  
  # Nest and test
  temp.df <- temp.df %>%
    nest() %>%
    mutate(aov = map(data, ~ aov(abundance ~ condition, data = .x)),
           summary = map(aov, tidy))
  
  # Unnest and FDR adjust
  temp.df <- temp.df %>%
    unnest(summary) %>%
    filter(term == "condition") %>%
    dplyr::rename(., P = p.value) %>%
    mutate(FDR = p.adjust(P, method = "BH", n = length(P)))
  
  # Join to data
  temp.df <- temp.df %>%
    select(1, P, FDR) %>%
    left_join(df, ., by = variable)
  
  # Exit
  return(temp.df)
  
}


calculate_fc <- function(data, group.compare){
  for (x in group.compare){
    temp.data <- data
    
    temp.fc <- rowMeans(temp.data[x[[2]]]) - rowMeans(temp.data[x[[1]]])

    # Add to dataframe
    data <- cbind(data, temp.name = temp.fc)
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- data %>%
      select(x[[1]][1]) %>%
      names() %>%
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
    temp2 <- data %>%
      select(x[[2]][1]) %>%
      names() %>%
      str_split(., "-", 2) %>% .[[1]] %>% .[1]
    
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
  
  # Rename Z-score columns
  names(temp.df)[-1] <- names(temp.df)[-1] %>% paste0(., "_Z")
  
  temp.df <- temp.df %>%
    select(-1) %>%
    
    mutate(Cluster = dist(.) %>%
             hclust() %>%
             cutree(k)) %>%
    
    mutate(Cluster = LETTERS[Cluster]) %>%
    
    bind_cols(temp.df[1], .)
  
  # Exit
  return(temp.df)
  
}
