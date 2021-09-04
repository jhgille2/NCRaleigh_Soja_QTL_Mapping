#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Phenotype_Data
read_clean <- function(Phenotype_Data) {

  # The sheets with the phenotype data
  Pheno_sheets <- excel_sheets(Phenotype_Data) %>% 
    str_subset("redo")
  
  # The location for each sheet of phenotype data
  Pheno_Locs <- map_chr(Pheno_sheets, function(x) str_split(x, " ", simplify = TRUE)[[1]])
  
  # Read in all the sheets
  AllPheno <- map(Pheno_sheets, function(x) read_excel(Phenotype_Data, sheet = x)) %>% 
    set_names(Pheno_Locs)
  
  # A function to clean each sheet
  Clean_sheet <- function(sheet, loc_name){
    
    sheet %<>% 
      clean_names() %>% 
      # Select a subset of the columns that will be used in the analysis
      select(sample_name, 
            # code, 
             wt_mg, 
             percent_nitrogen, 
             percent_carbon, 
             percent_sulfur) %>% 
      # Add the location name, and a plot number derived from the sample name
      mutate(loc = loc_name,
             plot = readr::parse_number(sample_name)) %>% 
      # Remove any sample names without a code (standards and controls in the spreadsheet)
      dplyr::filter(!is.na(code)) %>% 
      group_by(code) %>% 
      arrange(plot) %>% 
      # Add a rep number based on the plot numbers
      mutate(rep = 1:n()) %>% 
      ungroup() %>% 
      arrange(loc, code, rep) %>%
      mutate(n_s_percent_ratio = percent_nitrogen/percent_sulfur)
    
    return(sheet)
  }
  
  # Read in and clean each sheet, then combine into a single dataframe
  Combined_pheno <- map2(AllPheno, names(AllPheno), Clean_sheet) %>% 
    reduce(bind_rows) %>% 
    mutate(across(c(code, loc, rep), factor))
  
  return(Combined_pheno)
}
