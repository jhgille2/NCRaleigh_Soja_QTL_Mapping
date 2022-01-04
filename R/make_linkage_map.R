#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param genodata
#' @param phenodata
#' @param F.generation
#' @param missing_threshold
#' @param segdistortion_threshold
#' @return
#' @author 'Jay
#' @export
make_linkage_map <- function(genodata = Cleaned_Genotypes, phenodata =
                             Phenotype_Export, F.generation = 4,
                             missing_threshold = 0.05, segdistortion_threshold
                             = 0.001) {

  # Read in the cross
  snpData_cross <- read.cross(format = "csvsr", 
                              genfile = genodata, 
                              phefile = phenodata, 
                              F.gen = F.generation)
  
  # Remove colocated markers, markers with significant segregation distortion, 
  # and those with more than 5% missing data
  snpData_cross          <- pullCross(snpData_cross, type = "co.located")
  snpData_cross_segdist  <- pullCross(snpData_cross, type = "seg.distortion", pars = list(seg.thresh = 0.001))
  snpData_cross_noMiss_5 <- pullCross(snpData_cross_segdist, type = "missing", pars = list("miss.thresh" = 0.05))

  # Estimate the linkage map and keep the markers in their physical order
  snpData_mst_pass1_anchor   <- mstmap.cross(snpData_cross_noMiss_5, id = "id", anchor = TRUE)
  
  # Subset the map and keep only the markers that group to the main chromosomes
  lg_remainders <- as.numeric(chrnames(snpData_mst_pass1_anchor)) %% 1
  
  snpData_subsetted <- subsetCross(snpData_mst_pass1_anchor, chrnames(snpData_mst_pass1_anchor)[which(lg_remainders <= 0.11)])
  snpData_subsetted <- jittermap(snpData_subsetted)
  
  return(snpData_subsetted)
}
