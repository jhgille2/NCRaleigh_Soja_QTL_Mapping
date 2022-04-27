#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param cross
qtl_map_regression <- function(LinkageMap) {

  # Get the numeric indices for each phenotype
  get_pheno_indices <- function(pheno_name){
    
    pheno_indices <- str_detect(phenames(LinkageMap), pheno_name) %>% which()
    
    return(pheno_indices)
  }
  
  # Apply the function to get a list of numeric indices for each phenotype
  phenos <- c("n_s", "c_n", "sulfur", "nitrogen", "carbon")   
  all_pheno_indices <- map(phenos, get_pheno_indices) %>% 
    set_names(phenos)
  
  # A function to do both permutations and mapping for a set of phenotype
  # indices as returned by the function above
  perm_and_map <- function(pheno_indices){
    
    # zzzpermutations and mapping with marker regression for each phenotype
    perm <- scanone(LinkageMap, 
                    pheno.col = pheno_indices, 
                    method = "mr", 
                    n.perm = 1000, 
                    n.cluster = 10)
    
    scan <- scanone(LinkageMap, 
                    pheno.col = pheno_indices, 
                    method = "mr")
    
    # Return both as a list
    return(list(perms = perm, 
                result = scan))
  }
  
  # Apply this function to the lisst of indices to get a list of marker regression results
  all_mr_results <- map(all_pheno_indices, perm_and_map)

  # And return this list
  return(all_mr_results)
}
