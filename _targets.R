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
             clean_genotypes(genofile     = Genotype_data, 
                             new_genofile = New_Genotype_data,
                             parentA      = "2104", 
                             parentB      = "2105")), # Soja 
  
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
  # tar_target(CrossData,
  #            read_to_cross(genoFile  = Genotype_Export,
  #                          phenoFile = Phenotype_Export)),
  # 
  # # Remove genotypes with no phenotype data and convert the cross to a bcsft object with F.gen number of selfing generations
  # bcsft_mapped <- tar_map(
  #   
  #   unlist = FALSE,
  #   values = tibble(fgen = c(4, 5)), 
  #   
  #   tar_target(bcsftCross, 
  #              create_bcsft(Cross = CrossData, fgen))
  # ),
  # 
  # # Combine the bcsft crosses into a list
  # tar_combine(bcsft_combined,
  #             bcsft_mapped[[1]],
  #             command = list(!!!.x), 
  #             iteration = "list"),
  # 
  # # Some marker based quality control for the crosses, return a tibble of summary stats for 
  # # the cross after applying the filters
  # cross_qc_mapped <- tar_map(
  # 
  #   unlist = FALSE,
  #   values = expand_grid(miss = c(0.1, 0.05, 0.025),
  #                        seg.dist = c(0.1, 0.05, 0.025)),
  # 
  #   tar_target(bcsft_QC,
  #              apply_marker_filters(bcsft_combined, miss, seg.dist), 
  #              pattern = map(bcsft_combined))
  # 
  #   
  # ),
  # 
  # # Combine all the summaries into one table. 
  # tar_combine(all_snp_QC, 
  #             cross_qc_mapped[[1]], 
  #             command = dplyr::bind_rows(!!!.x)),
  # 
  # # Render the notebook and writeup
  # tar_render(Writeup, 
  #            here("docs", "Notebook.Rmd")),
  # 
  # tar_render(Journal, 
  #            here("docs", "Journal.Rmd"))

)
