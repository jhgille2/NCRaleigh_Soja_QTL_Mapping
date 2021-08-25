#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param filepath
#' @param genodata
export_genotype <- function(filepath = here("data", "rqtl_genotypes"), genodata
                            = ABH_Genotypes) {
  
  # The chromosomes for each marker
  MarkerChr <- as.numeric(str_sub(rownames(genodata), 3, 4))
  
  # The position of each marker
  MarkerPos <-  rownames(genodata) %>% 
    strsplit("_") %>% 
    map_chr(., function(x) x[[2]]) %>% 
    as.numeric()

  # A tibble to hold the name, chromosome, and position of each marker
  MarkerTable <- data.frame(mname = rownames(genodata),
                            chr = MarkerChr, 
                            pos = MarkerPos)
  
  # Bind this marker postion table to the genotype data 
  # (needed for the csvsr format)
  AllData <- cbind(MarkerTable, genodata)
  
  colnames(AllData) <- c(c("id", "", ""), colnames(genodata))
  
  # Export the data to a .csv file
  write.table(AllData, 
              file = filepath, 
              sep = ",", 
              row.names = FALSE, 
              col.names = TRUE)
    
  return(filepath)
}
