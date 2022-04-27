#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param oldgenodata
#' @param newgenodata
#' @return
#' @author 'Jay
#' @export
prep_genotype_data <- function(oldgenodata = Genotype_Data, newgenodata =
                               New_Genotype_data, 
                               sample_remeasures = Genotype_data_24_samples) {

  genodata1 <- vroom::vroom(oldgenodata, skip = 9)
  genodata2 <- vroom::vroom(newgenodata, skip = 9)
  genodata3 <- vroom::vroom(sample_remeasures, skip = 9)
  
  # Join the three sets of genotype data, excluding the samples that were remeasured 
  # from the first two sets of data before joining to the third
  snpData <- left_join(genodata1, genodata2, by = "...1") %>% 
    select(-any_of(colnames(genodata3)[2:ncol(genodata3)])) %>% 
    left_join(., genodata3, by = "...1") %>%
    column_to_rownames("...1") %>% 
    relocate(c("2104", "2105"))
  
  # Functions to get the chromosome and bp position from a snp name
  get_lg_from_snp_name <- function(snp_name)
  {
    snp_name %>% 
      str_split(., pattern = "_") %>% 
      pluck(1, 1) %>%
      parse_number() %>% 
      as.numeric()
  }
  
  get_bp_from_snp_name <- function(snp_name)
  {
    snp_name %>% 
      str_split(., pattern = "_") %>% 
      pluck(1, 2) %>%
      parse_number() %>% 
      as.numeric()
  }
  
  # Only keep numeric samples (the ones in this population)
  mp_samples <- sort(as.numeric(colnames(snpData))) %>% as.character()
  
  snpData <- snpData[, mp_samples] %>% 
    relocate(c("2104", "2105"))
  

  # Recode genotypes to ABH format according to the parental genotypes  
  hetMatrix <- snpData == "AB"
  snpData[snpData == snpData[, "2104"]] <- "A"
  snpData[hetMatrix] <- "H"
  snpData[snpData == "--"] <- "-"
  snpData[snpData == "BB"] <- "B"
  snpData[snpData == "AA"] <- "B"
  

  
  # Create chr and position columns following the chromosome assignment
  # and physical position of the SNPs
  snpData %<>% 
    mutate(chr = map_dbl(rownames(.), get_lg_from_snp_name), 
           pos = map_dbl(rownames(.), get_bp_from_snp_name)) %>% 
    relocate(chr, pos) %>% 
    rownames_to_column(var = "id")
  
  # Formatting to follow the r/qtl genotype file format
  colnames(snpData)[2:3] <- ""
  
  # Write the data to a csv file and return the path to this file
  write.table(snpData, file = here("data", "fullgenodata.csv"), row.names = FALSE, col.names = TRUE, sep = ",")
  return(here("data", "fullgenodata.csv"))
}
