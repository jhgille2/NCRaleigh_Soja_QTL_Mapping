#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Read_and_clean
make_eda_plots <- function(Read_and_clean) {

  cleanData_long <- Read_and_clean %>% 
    pivot_longer(cols = c(percent_nitrogen, percent_carbon, percent_sulfur, n_s_percent_ratio), names_to = "pheno")
  
  # A boxplot of phenotype values by location
  Pheno_boxplots <- ggplot(cleanData_long, aes(x = value, y = loc), fill = loc) + 
    geom_boxplot() +  
    theme_hc() + 
    facet_wrap(~pheno, scales = "free") + 
    coord_flip()
  
  # Histograms of phenotype values, faceted by phenotype and location

  Pheno_histograms <- ggplot(cleanData_long, aes(x = value)) + 
    geom_histogram() + 
    theme_hc() + 
    facet_grid(rows= vars(loc), cols = vars(pheno), scales = "free_x")
  
  # qq plots
  Pheno_qqplots <- ggplot(cleanData_long, aes(sample = value)) + 
    stat_qq() +
    stat_qq_line() + 
    theme_hc() + 
    facet_grid(rows= vars(loc), cols = vars(pheno), scales = "free")
  
  # Combine everything in a list and return
  allPlots <- list("boxplot"    = Pheno_boxplots, 
                   "histograms" = Pheno_histograms, 
                   "qqplots"    = Pheno_qqplots)
  
  return(allPlots)
}
