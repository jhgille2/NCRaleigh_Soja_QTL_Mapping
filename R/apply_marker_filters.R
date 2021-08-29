#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bcsft_mapped
#' @param miss
#' @param seg.dist
apply_marker_filters <- function(Cross, miss, seg.dist) {

  pull_coloc         <- pullCross(Cross, type = "co.located")
  pull_missing       <- pullCross(pull_coloc, type = "missing", pars = list(miss.thresh = miss))
  pull_segdistortion <- pullCross(pull_missing, type = "seg.distortion", pars = list(seg.thresh = seg.dist))
  
  # Return a tibble with some summary stats
  CrossSummary <- summary(pull_segdistortion)
  
  res <- tibble(misspct = miss, 
                seg     = seg.dist, 
                nmar    = sum(nmar(pull_segdistortion)), 
                F.gen   = CrossSummary$cross.scheme[[2]], 
                summary = list(CrossSummary))
  
  return(res)
}
