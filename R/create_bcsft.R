#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Cross
#' @param fgen
create_bcsft <- function(Cross = CrossData, fgen = 5) {

  # Remove individuals that don't have phenotype measurements
  PhenoSamples <- subset(Cross, ind = !is.na(pull.pheno(Cross)[, 2]))
  
  
  bcsftCross <- convert2bcsft(PhenoSamples, F.gen = fgen, map.function = "kosambi")
  
  mstmap_init <- mstmap.cross(bcsftCross, id = "id", p.value = 2)

  return(mstmap_init)
}
