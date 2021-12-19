n#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Genotype_data
clean_genotypes <- function(genofile = Genotype_data, new_genofile = New_Genotype_data, parentA = "2104", parentB = "2105") {
  
  # Read in the genotype matrices from the ols and the new data and then
  # join them together
  snpData     <- vroom::vroom(genofile, skip = 9)
  snpData_new <- vroom::vroom(new_genofile, skip = 9)
  
  snpData <- left_join(snpData, snpData_new, by = "...1") %>% 
    column_to_rownames("...1")

  ## Section: Helper functions for parent SNP QC
  ##################################################
  
  # A function to compare the genotypes of two sampels and return the indices where
  # they have the same genotypes
  compare_two_genotypes <- function(sample1_name = "2104", sample2_name = "2105", genotype_matrix = NULL)
  {
    # Get the genotypes for each sample
    sample1_genotypes <- genotype_matrix[, sample1_name]
    sample2_genotypes <- genotype_matrix[, sample2_name]
    
    # Get the indices of where the two samples have the same genotype
    same_genotype_indices <- which(sample1_genotypes == sample2_genotypes)
    
    return(same_genotype_indices)
  }
  
  # A function that returns the indices of a genotype that have missing genotypes
  find_missing_genotypes <- function(sample_name = "2104", genotype_matrix = NULL, missing_character = "--")
  {
    missing_indices <- which(genotype_matrix[, sample_name] == missing_character)
    
    return(missing_indices)
  }
  
  
  # A function to find the indices of heterozygous genotypes
  find_heterozygous_snps <- function(sample_name = "2104", genotype_matrix = NULL)
  {
    # A function that returns TRUE if a string has multiple characters
    is_het <- function(genotype){
      
      heterozygous <- sum(!!str_count(genotype, LETTERS)) > 1
      
      return(heterozygous)
      
    }
    
    # The indices of the heterozygous genotypes
    samplegeno       <- genotype_matrix[, sample_name]
    heterozygous_ids <- map_lgl(samplegeno, is_het) %>% which()
    
    return(heterozygous_ids)
  }
  
  # A function that combines the above functions and returns a vector of snp
  # indices to exclude
  snps_to_exclude <- function(sample1_name = "2104", sample2_name = "2015", genotype_matrix = NULL, missing_character = "--")
  {
    # The indices of the SNPs that are the same in each parent
    same_genotype_indices <- compare_two_genotypes(sample1_name = sample1_name, sample2_name = sample2_name, genotype_matrix = genotype_matrix)
    
    # The indices of the snps with missing data in both parents
    missing_snps_sample1 <- find_missing_genotypes(sample_name = sample1_name, genotype_matrix = genotype_matrix)
    missing_snps_sample2 <- find_missing_genotypes(sample_name = sample2_name, genotype_matrix = genotype_matrix)
    
    # The indices of the heterozygous SNPs in each sample
    het_snps_sample1 <- find_heterozygous_snps(sample_name = sample1_name, genotype_matrix = genotype_matrix)
    het_snps_sample2 <- find_heterozygous_snps(sample_name = sample2_name, genotype_matrix = genotype_matrix)
    
    # Find the unique SNPs from all of these conditions
    exclude_snps <- unique(c(same_genotype_indices, 
                             missing_snps_sample1, missing_snps_sample2, 
                             het_snps_sample1, het_snps_sample2))
    
    return(exclude_snps)
  }
  
  ##########################################################################
  
  
  # Combine all the indices to remove and then remove those SNPs
  all_remove_ids <- snps_to_exclude(sample1_name      = parentA, 
                                    sample2_name      = parentB, 
                                    genotype_matrix   = snpData, 
                                    missing_character = "--")
  
  # Remove these SNPs from the full data
  snpData <- snpData[-all_remove_ids, ]
  
  valid_heterozygotes <- function(){
    parentA_alleles <- map_chr(snpData[, parentA], function(x) str_sub(x, 1, 1))
    parentB_alleles <- map_chr(snpData[, parentB], function(x) str_sub(x, 1, 1))
    
    het1 <- paste0(parentA_alleles, parentB_alleles)
    het2 <- paste0(parentB_alleles, parentA_alleles)
    
    return(list(het1 = het1, het2 = het2))
  }
  
  allhets <- valid_heterozygotes()
  
  browser()
  
  # Convert genotypes to A, B, and H following the parental genotypes
  snpData[snpData == snpData[, parentA]] <- "A"
  snpData[snpData == snpData[, parentB]] <- "B"
  
  # Replace "--" with NA
  snpData[snpData == "--"] <- NA
  
  return(snpData)
}
