#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param genoFile
#' @param phenoFile
read_to_cross <- function(genoFile = Genotype_Export, phenoFile =
                          Phenotype_Export) {

  read.cross(format    = "csvsr", 
             genfile   = genoFile, 
             phefile   = phenoFile) # Check how appropriate it is to treat this as a selfed RI population

}
