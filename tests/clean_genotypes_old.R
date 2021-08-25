n#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Genotype_data
clean_genotypes <- function(genofile = Genotype_data, parentA = "1926", parentB = "2002", sojaFile = SojaParentGeno, snpConversion = SNP_Conversion) {

  # Read in the soja parent genotypes and the SNP name conversion table
  sojaGeno   <- read.vcfR(sojaFile)
  snpconvert <- vroom::vroom(snpConversion)
  
  
  # Read in the genotype matrix data
  snpData <- vroom::vroom(genofile, skip = 9) %>% 
    rename(secondaryIdentifier = "...1")
  
  # Extract the genotypes from the vcf object
  snpData <- extract.gt(sojaGeno, return.alleles = TRUE) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "SNP_ID") %>% 
    # remove slashes from genotypes
    mutate(PI424025B = str_remove(PI424025B, "\\/"), 
           PI424025B = replace_na(PI424025B, "--")) %>% 
    # Join with the snp conversion table to add the long snp names that are used in the chip
    left_join(snpconvert, by = c("SNP_ID" = "primaryIdentifier")) %>% 
    select(secondaryIdentifier, PI424025B) %>%
    # Replace the PI name with the code name that is used in the phenotype data
    rename(`2105` = PI424025B) %>%
    # Remove any SNPs that didn't have a name in the needed format
    filter(!is.na(secondaryIdentifier)) %>% 
    # Join with the chip data following the SNP names
    right_join(snpData, by = "secondaryIdentifier") %>%
    # Convert the SNP names bac to rownames of the data
    column_to_rownames(var = "secondaryIdentifier")
  
  # Some conditions to filter SNPs
  # 1. SNPs that have the same identity in each parent
  # 2. SNPs that are missing in either parent
  # 3. SNPs that are heterozygous in either parent
  
  # The genotype vectors for each parent
  parentA_geno <- snpData[, parentA]
  parentB_geno <- snpData[, parentB]
  
  # The indices of the SNPs with the same allelic identity in each parent
  same_allele_ids <- which(parentA_geno == parentB_geno)
  
  # The indices of the SNPs with missing data in either parent
  missing_ids <- which(parentA_geno == "--" | parentB_geno == "--")
  
  # A function that returns TRUE if a string has multiple characters
  is_het <- function(genotype){
    
    heterozygous <- sum(!!str_count(genotype, LETTERS)) > 1
    
    return(heterozygous)
    
  }
  
  het_parentA <- map_lgl(parentA_geno, is_het) %>% which()
  het_parentB <- map_lgl(parentB_geno, is_het) %>% which()
  
  heterozygous_ids <- unique(c(het_parentA, het_parentB))
  
  # Combine all the indices to remove and then remove those SNPs
  all_remove_ids <- unique(c(same_allele_ids, missing_ids, heterozygous_ids))
  
  snpData <- snpData[-all_remove_ids, ]
  
  # A function to get the heterozygous genotypes from parent genotypes
  valid_heterozygotes <- function(){
    parentA_alleles <- map_chr(snpData[, parentA], function(x) str_sub(x, 1, 1))
    parentB_alleles <- map_chr(snpData[, parentB], function(x) str_sub(x, 1, 1))
    
    het1 <- paste0(parentA_alleles, parentB_alleles)
    het2 <- paste0(parentB_alleles, parentA_alleles)
    
    return(list(het1 = het1, het2 = het2))
  }
  
  allhets <- valid_heterozygotes()
  
  snpData[snpData == snpData[, parentA]] <- "A"
  snpData[snpData == snpData[, parentB]] <- "B"
  
  # Find all heterozygous cells
  snpData[snpData == allhets$het1] <- "H"
  snpData[snpData == allhets$het2] <- "H"
  
  # Replace "--" with NA
  snpData[snpData == "--"] <- NA
  
  return(snpData)
}
