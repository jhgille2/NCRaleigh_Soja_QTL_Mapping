tar_load(LinkageMap)


all_perms <- scanone(cross = LinkageMap, 
                     pheno.col = 2:nphe(LinkageMap), 
                     method = "mr", 
                     n.perm = 1000, 
                     n.cluster = 10)

all_marker_regression <- scanone(cross = LinkageMap, 
                                 pheno.col = 2:nphe(LinkageMap), 
                                 method = "mr")


summary(all_marker_regression, 
        perms = all_perms, 
        format = "allpeaks")

# Indices for each phenotype
get_pheno_indices <- function(pheno_name){
  
pheno_indices <- str_detect(phenames(LinkageMap), pheno_name) %>% which()

return(pheno_indices)
}

phenos <- c("n_s", "sulfur", "nitrogen", "carbon")   
all_pheno_indices <- map(phenos, get_pheno_indices) %>% 
  set_names(phenos)

# A function to do both permutations and mapping for a set of phenotypes
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

all_mr_results <- map(all_pheno_indices, perm_and_map)
