##################################################
## Project: Raleigh x Soja QTL Mapping
## Script purpose: Helper functions
## Date: 2022-02-07
## Author: Jay Gillenwater
##################################################

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
