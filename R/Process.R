reduce_features <- function(df){
  temp.df <- df %>%
    # Summarize features with identical peptides - different accessions
    group_by(`#`, Sequence, Modifications, Score) %>%
    top_n(1, `Unique peptides`) %>%
    top_n(1, `Confidence score`) %>%
    dplyr::slice(1) %>%
    
    # Summarize features with multiple peptides
    group_by(`#`) %>%
    top_n(1, Score) %>%
    dplyr::slice(1) %>%
    ungroup()
  
}


get_identifier <- function(df, database, mod = "Phospho"){
  # Separate to one modification per row
  temp.df <- df %>%
    select(`#`, Accession, Sequence, Modifications) %>%
    separate_rows(Modifications, sep = "\\|") %>%
    filter(str_detect(Modifications, mod))
  
  # Get the modification position and residue from the peptide
  temp.df <- temp.df %>%
    mutate(position = str_extract(Modifications, "(?<=\\[)(.*)(?=\\])") %>%
             as.numeric()) %>%
    mutate(residue = str_sub(Sequence, start = position, end = position))
  
  # Get the modification site from the protein
  temp.df <- temp.df %>%
    mutate(locate = str_locate(database[Accession], Sequence)[, 1]) %>%
    mutate(site = locate + position - 1)
  
  # Group by feature and build identifier
  temp.df <- temp.df %>%
    group_by(`#`) %>%
    mutate(Identifier = paste(residue, site, sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--")) %>%
    ungroup()
  
  # Join identifier column to input data
  temp.df %>%
    select(`#`, Identifier) %>%
    distinct() %>%
    inner_join(df, ., by = "#")
  
}


reduce_identifiers <- function(df, samples){
  # Reduce to unique identifiers represented by highest score
  temp.df <- df %>%
    group_by(Identifier) %>%
    top_n(1, Score) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Column-wise sum
  temp.df2 <- df %>%
    select(Identifier, samples) %>%
    gather(sample, abundance, -1) %>%
    mutate(sample = factor(sample, levels = df[, samples] %>% names())) %>%
    
    group_by(Identifier, sample) %>%
    summarize(sum = sum(abundance)) %>%
    
    spread(sample, sum) %>%
    ungroup()
  
  # Replace summed abundance in ordered dataset
  temp.df[, samples] <- temp.df2[, -1]
  
  return(temp.df)

}


filter_redox <- function(df, reduced = "Nethylmaleimide"){
  temp.df <- df %>%
    mutate(Modifications = replace_na(Modifications, "")) %>%
    filter(., str_count(Sequence, "C") > str_count(Modifications, reduced))
  
}


get_identifier_redox <- function(df, database, reduced = "Nethylmaleimide"){
  # Map integer position of modification
  temp.df <- df %>%
    mutate(Modified = str_split(as.character(Modifications), "\\|") %>%
             lapply(., str_subset, reduced) %>%
             lapply(., str_extract, "(?<=\\[)(.*)(?=\\])") %>%
             lapply(., as.numeric))
  
  # Map position of residue to peptide
  temp.df <- temp.df %>%
    rowwise() %>%
    mutate(Residue = str_locate_all(Sequence, "C")[[1]][, 1] %>% list())
  
  # Map unmodified position of residue on peptide
  temp.df <- temp.df %>%
    mutate(Residue = setdiff(Residue, Modified) %>% list())
  
  # Map position of peptide to protein
  temp.df <- temp.df %>%
    mutate(Start = str_locate(database[Accession], Sequence)[[1]])
  
  # Build identifier
  temp.df <- temp.df %>%
    mutate(Identifier = (Residue + Start - 1) %>%
             paste("C", ., sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--"))
  
  # Return clean dataframe
  temp.df <- temp.df %>%
    select(-Modified, -Residue, -Start) %>%
    ungroup()

}
  
  
get_identifier_redox2 <- function(df, database, reduced = "Nethylmaleimide"){
  # TODO
  
  # Locate each Cys residue and separate into rows
  temp.df <- df %>%
    select(`#`, Accession, Sequence, Modifications) %>%
    rowwise() %>%
    mutate(Cys = str_locate_all(Sequence, "C")[[1]][, 1] %>% paste(., collapse = "-")) %>%
    separate_rows(Cys, sep = "-")
  
  # Separate modifications and remove blocked Cys
  temp.df %>%
    separate_rows(Modifications, sep = "\\|")
  
}


simulate_data <- function(data = 1000, p = 100, fold = c(2, 4), num.rep = 3, cv = 20){
  # Initiate RNG state
  set.seed(123)
  
  temp.df <- tibble(Accession = c(str_c("N", c(1:(data - p)), sep = "_"),
                                  str_c("P", c(1:p), sep = "_")))
  
  temp.df <- temp.df %>%
    mutate(Range = rnorm(n = n(), mean = 8, sd = 3))
  
  fold2 <- c(1, fold)
  
  for (x in 1:length(fold2)){
    for (y in 1:num.rep){
      temp.name <- LETTERS[x] %>%
        paste(., y, sep = "-")
      
      temp.df <- temp.df %>%
        mutate(!!temp.name := if_else(str_detect(Accession, "P_"),
                                      rnorm(n = Range, mean = Range + log2(fold2[x]), sd = abs(Range * (cv / 100))),
                                      rnorm(n = Range, mean = Range, sd = abs(Range * (cv / 100)))))
    }
  }
  temp.df <- temp.df %>%
    select(-Range)
  
  temp.df %>%
    mutate(sd = abs(Range * (cv / 100)))
  
}