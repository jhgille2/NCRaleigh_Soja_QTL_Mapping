#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Marker_regression_results
summarise_makrer_regression <- function(Marker_regression_results) {

  # A function to get a table summarizing the significant qtl for a 
  # given alpha value for an element of the Marker_regression_results list
  summarise_sig_qtl <- function(mr_result, alpha = 0.1){
    
    # First, use the permutations to declare significant qtl
    by_loc_summary <- summary(mr_result$result, 
                              format = "tabByCol", 
                              perms = mr_result$perms, 
                              alpha = alpha)
    
    # Only keep the locations which had a significant QTL
    has_sig_qtl <- map_dbl(by_loc_summary, nrow) > 0
    by_loc_summary <- by_loc_summary[has_sig_qtl]
    
    # by_loc_summary is a named list where each list element
    # had the significant qtl for the environemnt specified by the
    # name of the list element. I want to add a location column to 
    # each result dataframe to hold onto this location
    add_location <- function(summary_table, loc_name){
      
      sig_markers <- rownames(summary_table)
      
      summary_table %<>% 
        as_tibble() %>%
        mutate(Environment = loc_name, 
               `Linked Marker` = sig_markers)
      
      return(summary_table)
    }
    
    # Now, apply this function to the named list of sumnmary tables
    # and then merge all the summary tables
    by_loc_summary %<>%
      map2(., names(.), add_location) %>% 
      unname() %>%
      reduce(bind_rows)
    
    return(by_loc_summary)
  }
  
  # A function to parse the "location" column. This column has the 
  # long and cumbersome names that I used to keep track of the location
  # and the trait together. Now, simplify this into a more readible 
  # format
  extract_location_name <- function(longLocName){
    shortLocName <- ifelse(str_detect(longLocName, "LA_"), "LA", 
                           ifelse(str_detect(longLocName, "SAN_"), "SAN", 
                                  ifelse(str_detect(longLocName, "MS_"), "MS", 
                                         ifelse(str_detect(longLocName, "CLA_"), "CLA", 
                                                ifelse(str_detect(longLocName, "SC_"), "SC", "Average")))))
    
    return(shortLocName)
  }
  
  # The mr results list is named by phenotype. Add a column to the 
  # summary table that keeps this phenotype name
  add_phenotype_name <- function(summary_table, pheno_name){
    summary_table %<>% 
      mutate(Phenotype = pheno_name)
    
    return(summary_table)
  }
  
  # Now, apply this function to all the list elements and bind everything to
  # one big dataframe
  Marker_summary_table <- map(Marker_regression_results, summarise_sig_qtl) %>% 
    map2(., names(.), add_phenotype_name) %>% 
    reduce(bind_rows) %>% 
    mutate(Environment = map_chr(Environment, extract_location_name)) %>% 
    select(Phenotype, Environment, chr, pos, `Linked Marker`, ci.low, ci.high, lod) %>% 
    arrange(Phenotype, Environment, chr, pos)
  
  return(Marker_summary_table)
}
