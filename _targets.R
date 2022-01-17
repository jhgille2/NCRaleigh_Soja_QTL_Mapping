## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

# Parallelism with future
plan(callr)

tar_option_set(error = "workspace")

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
             here("data", "snp_matrix_AB.csv"), 
             format = "file"),
  
  # A SNP file for some samples that were genotyped later
  tar_target(New_Genotype_data, 
             here("data", "new_samples_AB.csv"), 
             format = "file"),
  
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
  tar_target(Cleaned_Genotypes, 
             prep_genotype_data(oldgenodata = Genotype_data, 
                                newgenodata = New_Genotype_data), 
             format = "file"),
  
  ## Section: Export files for r/qtl
  ##################################################

  # Export the phenotype data to a file
  tar_target(Phenotype_Export,
             export_phenotype(filepath  = here("data", "rqtl_phenotypes.csv"),
                              phenodata = GenoMeans),
             format = "file"),
  
  # ## Section: Linkage mapping
  # ##################################################
  tar_target(LinkageMap, 
             make_linkage_map(genodata                = Cleaned_Genotypes, 
                              phenodata               = Phenotype_Export, 
                              F.generation            = 4, 
                              missing_threshold       = 0.05, 
                              segdistortion_threshold = 0.001)),
  
  ## Section: Document rendering
  ##################################################
  tar_render(progress_report, "doc/progress_report.Rmd")

)
