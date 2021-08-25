#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param Read_and_clean
remove_bad_samples <- function(Read_and_clean) {

  # Right now, it looks like only one sample should 
  # definitely be removed (LA PLOT: 255). This 
  # sample had a weight far higher than every other sample
  # and percent N, C, S far lower than every other sample. 
  # It seems likely that this is due to measurement error
  # rather than a legitimate sample as the other replicate
  # has values similar to those of the other entries. 
  #
  # This function will of course have to be adjusted if other
  # criteria are used in the future which identify additional
  # outliers. For now, it looks like this sample may be the only
  # (easily identifiable) outlier
  
  Cleaned_noOutlier <- Read_and_clean %>% 
    filter(sample_name != "LA PLOT: 255")

  return(Cleaned_noOutlier)
}
