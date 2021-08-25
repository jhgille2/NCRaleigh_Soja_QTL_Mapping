## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
  
  ## Section: Input files
  ##################################################

  # The phenotype data
  tar_target(Phenotype_Data, 
             here("data", "Pheno_five_locations.xlsx"), 
             format = "file"), 
  
  # The SNP file
  tar_target(Genotype_data, 
             here("data", "snp_matrix.csv"), 
             format = "file"),
  
  # # The file that has the genotypes for PI424025B (code 2105) from soybase
  # tar_target(SojaParentGeno, 
  #            here("data", "sojaParent.vcf"), 
  #            format = "file"), 
  # 
  # # A table that to convert between SNP ids
  # tar_target(SNP_Conversion, 
  #            here("data", "SNP_ConversionTable.csv"), 
  #            format = "file"),
  
  ## Section: Phenotype processing/EDA
  ##################################################
  
  # Read in and clean the phenotype data
  tar_target(Read_and_clean, 
             read_clean(Phenotype_Data)), 
  
  # Summary graphics
  tar_target(EDA_Plots, 
             make_eda_plots(Read_and_clean)), 
  
  # Summary tables
  tar_target(EDA_tables, 
             make_eda_tables(Read_and_clean)),
  
  # Remove a bad sample
  tar_target(RemoveErrors, 
             remove_bad_samples(Read_and_clean)),
  
  # Calculate the average and genotype LSMean 
  tar_target(GenoMeans, 
             calc_geno_means(RemoveErrors)),
  
  ## Section: Cleaning genotype data
  ##################################################
  
  # convert the raw genotypes to ABH format using the parent genotypes
  tar_target(ABH_Genotypes, 
             clean_genotypes(genofile = Genotype_data, 
                             parentA = "2104", 
                             parentB = "2105")), # Soja 
  
  ## Section: Export files for r/qtl
  ##################################################
  
  # Export genotype and phenotype data to seperate files in the r/qtl csvsr format
  tar_target(Genotype_Export,
             export_genotype(filepath = here("data", "rqtl_genotypes.csv"),
                             genodata = ABH_Genotypes),
             format = "file"),

  tar_target(Phenotype_Export,
             export_phenotype(filepath  = here("data", "rqtl_phenotypes.csv"),
                              phenodata = GenoMeans),
             format = "file"),

  # ## Section: Import to r/qtl cross
  # ##################################################
  tar_target(CrossData,
             read_to_cross(genoFile  = Genotype_Export,
                           phenoFile = Phenotype_Export)),
  
  # Render the notebook and writeup
  tar_render(Writeup, 
             here("docs", "Notebook.Rmd")),
  
  tar_render(Journal, 
             here("docs", "Journal.Rmd"))

)
