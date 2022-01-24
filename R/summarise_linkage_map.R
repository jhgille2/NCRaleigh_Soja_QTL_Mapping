#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param map
summarise_linkage_map <- function(lmap = LinkageMap) {

  # Make a tibble to hold information about the linkage map: 
  # For each chromosome/LG...
  #
  # 1. The number of markers
  # 2. The average marker spacing
  # 3. The overall LG size
  #
  
  # A function to get the average marker spacing for a linkage group
  average_marker_spacing <- function(lg)
  {
    avg_spacing <- lg %>% 
      diff() %>% 
      mean() %>% 
      round(2)
    
    return(avg_spacing)
  }
  
  # Make the summary table
  map_summary_table <- tibble(Chromosome                    = names(lmap$geno), 
                              `Number of markers`           = nmar(lmap), 
                              `Average marker spacing (cM)` = map_dbl(pull.map(lmap), average_marker_spacing), 
                              `LG size (cM)`                = c(chrlen(lmap)) %>% round(2))
  
  # Add a row that has the average of the numeric columns
  map_summary_table_withAverages <- bind_rows(map_summary_table, 
                                              map_dfc(map_summary_table, function(x) ifelse(is.numeric(x), mean(x), "Average")))
  
  return(map_summary_table_withAverages)
}
