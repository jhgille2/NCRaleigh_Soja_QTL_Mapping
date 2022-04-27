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
  
  # Genotype clones (No genetic clones in current data set)
  gc <- genClones(snpData_cross_noMiss_5, id = "id")
  
  # Remove clones
  # snpData_cross_noMiss_5 <- fixClones(snpData_cross_noMiss_5, 
  #                                     gc$cgd, 
  #                                     consensus = TRUE, 
  #                                     id = "id")
  
  # Remove any individuals with more than 10% missing marker data
  sg_miss <- statGen(snpData_cross_noMiss_5, id = "id", stat.type = "miss", bychr = FALSE)
  snpData_cross_noMiss_5 <- subset(snpData_cross_noMiss_5, ind = sg_miss$miss < sum(nmar(snpData_cross_noMiss_5)) * 0.1)

  # Estimate the linkage map and keep the markers in their physical order
  snpData_mst_pass1_anchor   <- mstmap.cross(snpData_cross_noMiss_5, 
                                             id = "id", 
                                             anchor = TRUE)
  
  # Chromosomes to keep 
  keep_chrs <- c("1", 
                 "2.1", 
                 "3", 
                 "4", 
                 "5", 
                 "6", 
                 "7.1", 
                 "8", 
                 "9.1", 
                 "10", 
                 "11", 
                 "12", 
                 "13.1", 
                 "14", 
                 "15", 
                 "16.1", 
                 "17.1", 
                 "18", 
                 "19", 
                 "20")
  
  # Subset the map and keep only the markers that group to the main chromosomes
  lg_remainders <- as.numeric(chrnames(snpData_mst_pass1_anchor)) %% 1
  
  snpData_subsetted             <- subsetCross(snpData_mst_pass1_anchor, keep_chrs)
  names(snpData_subsetted$geno) <- as.character(c(1:20))
  
  snpData_subsetted <- jittermap(snpData_subsetted)
  
  return(snpData_subsetted)
}
