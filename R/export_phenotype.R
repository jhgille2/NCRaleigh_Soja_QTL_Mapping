#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param filepath
#' @param phenodata
export_phenotype <- function(filepath = here("data", "rqtl_phenotypes"),
                             phenodata = GenoMeans) {
  
  # Get the mean data from the phenodata calculations list
  Phenotypes <- phenodata$average
  
  # Transpose this data to match the csvsr format for r/qtl
  Phenotypes %<>% 
    rename(id = code) %>% 
    t()
  
  # Write to a .csv file
  write.table(Phenotypes, 
              file      = filepath, 
              sep       = ",", 
              col.names = FALSE, 
              row.names = TRUE)

  # Return the path the file was saved to
  return(filepath)
}
