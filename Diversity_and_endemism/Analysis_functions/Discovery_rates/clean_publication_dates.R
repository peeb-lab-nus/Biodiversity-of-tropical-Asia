clean_publication_dates <- function(x, max = TRUE, first = FALSE){
  # Cleaning publication dates from the WCVP data and adapted to other taxa
  #
  # Arguments:
  #   x, chr, vector of publication dates (`first_published`)  
  #   max, logical, whether to choose the oldest year in the string (if FALSE, then minimum is taken)
  #   first, logical, whether to choose the first number in the string (overrides max = TRUE)
  #
  # Returns:
  #   numeric
  
  if(max == TRUE & first == TRUE) {
    warning("Both 'max' and 'first' == TRUE. Taking first number")
  }
  
  res <- vector()
  
  for(i in 1:length(x)){
    if(x[i] == "" | is.na(x[i])){
      res[i] <- NA
    } else {
      match_positions <- gregexpr("[0-9]+", text = x[i])
      numbers <- regmatches(x[i], match_positions)  
      numbers_vector <- as.numeric(unlist(numbers))
      if(first == FALSE) {
        if(max == TRUE){
          res[i] <- max(numbers_vector)
        } else {
          res[i] <- min(numbers_vector)
        }
      } else {
        res[i] <- numbers_vector[1]
      }
    }
  }
  return(res)
}