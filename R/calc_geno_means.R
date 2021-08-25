#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param RemoveErrors
calc_geno_means <- function(RemoveErrors) {

  # Phenotype averages by code
  code_averages <- RemoveErrors %>% 
    pivot_longer(cols = tidyselect::contains("percent")) %>% 
    group_by(code, name) %>% 
    summarise(pheno_avg = mean(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "name", values_from = pheno_avg)
  
  # Genotype marginal means for each phenotype
  
  # A function that fits a linear model and gets the marginal means from this model
  pheno_emmeans <- function(phenoData){
    
    marmeans <- lm(value ~ code*loc + rep:loc, data = phenoData) %>% 
      emmeans("code") %>% 
      as_tibble() %>% 
      select(code, emmean)
    
    return(marmeans)
  }
  
  pheno_emmeans <- RemoveErrors %>% 
    pivot_longer(cols = tidyselect::contains("percent")) %>%
    group_by(name) %>% 
    nest() %>% 
    mutate(mar_means = map(data, pheno_emmeans) )%>% 
    unnest(mar_means) %>% select(-one_of("data")) %>% 
    pivot_wider(names_from = "name", values_from = "emmean", names_prefix = "emmean_")
  
  # combine averages and marginal means in a single dataframs
  AllMeans <- left_join(code_averages, pheno_emmeans, by = "code")
  
  # Combine everything into a named list and return
  res <- list("average"  = code_averages, 
              "emmeans"  = pheno_emmeans, 
              "allmeans" = AllMeans)

  return(res)
}
