add_missingness <- function(df, raw, group){
  variable <- df %>%
    select(1) %>%
    names()
  
  # Filter raw data for StatLFQ passing observations
  temp.raw <- raw %>%
    semi_join(., df, by = variable)
  
  # Match column names
  names(temp.raw) <- names(df)[c(1, group %>% unlist())]
  
  # Convert to long-format and summarize for complete missingness
  temp.raw <- temp.raw %>%
    gather(sample, value, -1) %>%
    separate(sample,
             into = c("condition", "replicate"),
             sep = "-") %>%
    #select(-replicate) %>%
    group_by(!!as.name(variable), condition) %>%
    summarize(sum = sum(value), Missing = sum == 0) %>%
    filter(Missing == TRUE)
  
  # Add column with groups completely missing
  temp.raw <- temp.raw %>%
    group_by(!!as.name(variable)) %>%
    mutate(Missing = paste(condition, collapse = "-")) %>%
    distinct(., !!as.name(variable), .keep_all = TRUE) %>%
    select(variable, Missing)
  
  # Join missingness column onto data
  temp.data <- df %>%
    left_join(., temp.raw, by = variable)
  
  # Exit
  return(temp.data)
  
}


add_accession <- function(df){
  temp.df <- df %>%
    separate(Identifier,
             into = c("Accession", "Sites"),
             sep = "--",
             remove = FALSE) %>%
    select(-Accession, -Sites, everything())
  
}


add_uniprot <- function(df, path){
  # Parse "Accession" column for UniProt entry names
  temp.df <- df %>%
    mutate(Accession = str_split(Accession, pattern = "\\|", simplify = TRUE)[, 2])
    
  # Load UniProt annotation file
  annotation <- read_delim(path, delim = "\t", col_types = cols()) %>%
    dplyr::rename(Accession = Entry)
  
  # Join onto input data
  temp.df <- annotation %>%
    left_join(temp.df, ., by = "Accession")
  
}


pull_uniprot <- function(output = "Cr_uniprot_20190130_annotation.tsv"){
  path <- "https://www.uniprot.org/uniprot/?query=proteome%3AUP000006906&columns=id%2Centry%20name%2Cprotein%20names%2Cgenes%2Clength%2Cfeature(ACTIVE%20SITE)%2Cfeature(BINDING%20SITE)%2Ccomment(CATALYTIC%20ACTIVITY)%2Cfeature(METAL%20BINDING)%2Cgo(biological%20process)%2Cgo(cellular%20component)%2Cgo(molecular%20function)%2Ccomment(SUBCELLULAR%20LOCATION)%2Cfeature(DISULFIDE%20BOND)%2Ccomment(POST-TRANSLATIONAL%20MODIFICATION)%2Cdatabase(STRING)%2Cdatabase(KEGG)%2Cdatabase(GeneID)%2Cdatabase(KO)%2Cdatabase(Pfam)%2Ccomment(PATHWAY)%2Cdatabase(InterPro)%2Cdatabase(PhosphoSitePlus)&format=tab"
  
  # Pull table from UniProt URL
  temp.data <- read_delim(path, delim = "\t", col_types = cols())
  
  # Write table to local file
  #write_tsv(temp.data, output)
  
}
