filter_score <- function(df, threshold = 13){
  temp.df <- df %>%
    filter(Score != "---") %>%
    mutate(Score = Score %>% as.numeric()) %>%
    filter(Score > threshold)
  
}


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

remove_crap <- function(df){
  
crap_accessions <- c("P02769", "P00766", "P00767", "P00711", "Q7M135", "P00792", "P00791", "Q10735", "P30879",
                     "P04188", "P00760", "Q29463", "P00761", "P02662", "P02663", "P02666", "P02668", "P42212", 
                     "P81054", "P04745", "P13645", "P35527", "P04264", "P35908", "Q15323", "Q14532", "O76011", 
                     "Q92764", "O76013", "O76014", "O76015", "O76009", "Q14525", "Q14533", "Q9NSB4", "P78385", 
                     "Q9NSB2", "P78386", "O43790", "O77727", "P02534", "P25690", "P02539", "P15241", "P25691", 
                     "P02444", "P02445", "P02443", "P02441", "Q02958", "P02438", "P02439", "P02440", "P08131", 
                     "P26372", "O82803", "P15252", "P00004", "P00921", "P00330", "P00883", "P00698", "P68082", 
                     "P01012", "P00722", "P00366", "P02768", "P01008", "P08758", "P61769", "P55957", "P00915", 
                     "P00918", "P04040", "P07339", "P08311", "P01031", "P02741", "P00167", "P99999", "P01133", 
                     "P05413", "P06396", "P08263", "P09211", "P69905", "P68871", "P01344", "P10145", "P06732", 
                     "P00709", "P41159", "P61626", "P02144", "Q15843", "P15559", "P16083", "P01127", "P62937", 
                     "Q06830", "P01112", "P02753", "P00441", "P63165", "P12081", "P10636", "P10599", "P01375", 
                     "P02787", "P02788", "P51965", "O00762", "P63279", "P62979", "P32503")

temp.df <- df %>% 
  filter(!Protein.Group %in% crap_accessions)

}

database_parsing <- function(database){
  
as.data.frame(database)-> database

database2 <- database %>% 
  t() %>% 
  as.data.frame()

database2 <- database2 %>%
  rownames_to_column(var = "Accession") %>% 
  separate(Accession, 
           into = c("sp", "Accession"),
           by = "\\tr.") %>%  # ??? change the separation by = "" based on the leading characters before the database accessions
  select(1:2, "ProteinSequence" = V1)

}

get_identifier_FragPipe <- function(df, database, mod = "\\(UniMod\\:21\\)"){

# Get modification position on peptide #
  
  for (i in 1:nrow(df)) {
    current_sequence <- df$Modified.Sequence[i]
    all_positions <- c()
    
    # A while loop to handle sequences with multiple modifications
    while (str_count(current_sequence, mod) > 0) {
      # Find the starting position of the first modification
      position <- str_locate(current_sequence, mod) 
        
        # Store the position and add to the list
        all_positions <- c(all_positions, position[1])
        
        # Remove the first modification found
        current_sequence <- str_replace(current_sequence, mod, "")
    }
    
    # Paste all found positions into a single string
    df$positions[i] <- paste(all_positions, collapse = ";")
  }  
  
  # Separate to one modification per row
  temp.df <- df %>%
  select(Feature, Accession, Modified.Sequence, Sequence = Stripped.Sequence, positions) %>%
  separate_rows(positions, sep = "\\;")

  temp.df$positions <- as.numeric(temp.df$positions)
  

  
  temp.df <- temp.df %>% 
    mutate(positions = positions - 1) %>%
    mutate(residue = str_sub(Sequence, start = positions, end = positions))

# Get the modification site from the protein
temp.df <- temp.df %>%
  left_join(., database, by = "Accession") %>% 
  mutate(locate = str_locate(ProteinSequence, Sequence)[,1]) %>%
  mutate(site = locate + positions - 1)

# Group by feature and build identifier
temp.df <- temp.df %>%
  group_by(`Feature`) %>%
  mutate(Identifier = paste(residue, site, sep = "") %>%
           paste(., collapse = "-") %>%
           paste(Accession, ., sep = "--")) %>%
  ungroup()

# Join identifier column to input data
temp.df %>%
  select(`Feature`, Identifier, residue, site) %>%
  distinct() %>%
  inner_join(df, ., by = "Feature")

}

reduce_identifiers_FragPipe <-function(df, samples){
  # Reduce to unique identifiers represented by highest score
  temp.df <- df %>%
    group_by(Identifier) %>%
    top_n(1, Feature) %>%
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

get_phospho_localization <- function(df, database, Threshold = 0.75){
  
  temp.df <- df %>% 
    left_join(., pepm3, by = "Identifier") %>% 
    #select(Feature, Identifier, Accession, Sequence = Stripped.Sequence, localization =`STY:79.96633`, positions) %>% 
    mutate(Location = sapply(str_locate_all(Stripped.Sequence, "[SYT]"), function(x) paste(x[,"start"], collapse = ","))) %>% 
    mutate(Scores = sapply(str_extract_all(`STY:79.96633`, "(?<=\\()[0-9.]+(?=\\))"), paste, collapse = ",")) %>% 
    separate_rows(Location, Scores, sep = "\\,") %>% 
    separate_rows(positions, sep = "\\;")
  
  temp.df$positions <- as.numeric(temp.df$positions)
  temp.df$Scores <- as.numeric(temp.df$Scores)
  temp.df$Location <- as.numeric(temp.df$Location)
  
  
  temp.df <- temp.df %>% 
    mutate(positions = positions - 1) %>% 
    filter(Scores > Threshold) %>% 
    mutate(Scores = Scores*100) %>% 
    filter(positions == Location) %>% 
    mutate(residue = str_sub(Stripped.Sequence, start = positions, end = positions)) 
  
  
  temp.df <- temp.df %>%
    left_join(., database, by = "Accession") %>% 
    mutate(locate = str_locate(ProteinSequence, Stripped.Sequence)[,1]) %>%
    mutate(site = locate + positions - 1)
  
  # Group by feature and build identifier
  temp.df <- temp.df %>%
    group_by(`Feature`) %>%
    mutate(`Localization Score`= paste(residue, site, sep ="") %>% 
             paste(": ", Scores, "%", sep = "") %>% 
             paste(., collapse = ", ")) %>%
    ungroup()
  
# Join to input dataframe
  temp.df <- temp.df %>%
    select(Identifier, `Localization Score`) %>%
    distinct() %>% 
    inner_join(., df, by = "Identifier")
  
}
