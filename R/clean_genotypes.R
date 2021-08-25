n#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Genotype_data
clean_genotypes <- function(genofile = Genotype_data, parentA = "2104", parentB = "2105") {
  
  # The soja parent (2105) was not genotypes in this set. Instead, I'll recover it's genotype
  # from the NC-Raleigh genotype (2104) and infer it's genotype at each marker from the 
  # alternate allele at each marker. 
  #
  # Unfortunately this does mean that some quality control filters can't be done right now
  # like removing SNPs with missing data or were heterozygous in the soja parent.
  
  # Read in the genotype matrix data
  snpData <- vroom::vroom(genofile, skip = 9) %>% 
    column_to_rownames("...1")

  # The genotype vector for parent A
  parentA_geno <- snpData[, parentA]
  
  # The indices of the SNPs with missing data in this parent
  missing_ids <- which(parentA_geno == "--")
  
  # A function that returns TRUE if a string has multiple characters
  is_het <- function(genotype){
    
    heterozygous <- sum(!!str_count(genotype, LETTERS)) > 1
    
    return(heterozygous)
    
  }
  
  # The indices of the heterozygous genotypes in parent A
  heterozygous_ids <- map_lgl(parentA_geno, is_het) %>% which()
  
  # A function to get the alleles observed at a SNP
  GetSnpAlleles <- function(SNP){
    
    UniqueChars <- SNP %>% 
      unlist() %>% 
      strsplit("") %>% 
      unlist() %>% 
      unique()
    
    UniqueChars[UniqueChars != "-"]
  }
  
  # The alleles present at each marker
  AllAlleles <- apply(snpData, 1, GetSnpAlleles)
  
  # Monomorphic SNPS
  mono_snps <- which(map(AllAlleles, length) == 1)
  
  # Combine all the indices to remove and then remove those SNPs
  all_remove_ids <- unique(c(missing_ids, heterozygous_ids, mono_snps))
  
  snpData    <- snpData[-all_remove_ids, ]
  AllAlleles <- AllAlleles[-all_remove_ids]
  
  # The alleles for parent A
  parentA_Allele <- map_chr(snpData[, parentA], GetSnpAlleles)
  
  # A function that returns the allele at a locus that is not the parent A allele
  get_alt_allele <- function(parent_allele, all_alleles){
    return(all_alleles[all_alleles != parent_allele])
  }
  
  # The parent b haplotype
  parentB_haplo <- map2_chr(parentA_Allele, AllAlleles, get_alt_allele)
  
  # Double this haplotype to make the parent B genotype
  makeGeno <- function(haplo){
    return(paste0(haplo, haplo))
  }
  
  # The parent B genotype
  parentB_geno       <- map_chr(parentB_haplo, makeGeno)
  snpData[, parentB] <- parentB_geno
  
  # A function to get the heterozygous genotypes from parent genotypes
  valid_heterozygotes <- function(B_geno = parentB_geno){
    parentA_alleles <- map_chr(snpData[, parentA], function(x) str_sub(x, 1, 1))
    parentB_alleles <- map_chr(snpData[, parentB], function(x) str_sub(x, 1, 1))
    
    het1 <- paste0(parentA_alleles, parentB_alleles)
    het2 <- paste0(parentB_alleles, parentA_alleles)
    
    return(list(het1 = het1, het2 = het2))
  }
  
  allhets <- valid_heterozygotes()
  
  # Convert genotypes to A, B, and H following the parental genotypes
  snpData[snpData == snpData[, parentA]] <- "A"
  snpData[snpData == snpData[, parentB]] <- "B"
  
  # Find all heterozygous cells
  snpData[snpData == allhets$het1] <- "H"
  snpData[snpData == allhets$het2] <- "H"
  
  # Replace "--" with NA
  snpData[snpData == "--"] <- NA
  
  return(snpData)
}
