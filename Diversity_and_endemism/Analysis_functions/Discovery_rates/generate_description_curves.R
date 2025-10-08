
# ==================================================================================================
## DEFINE SPECIES ACCUMULATION CURVE FUNCTION
# ==================================================================================================

generate_description_curves <- function(x,
                                        ladderise = TRUE,
                                        min_year = NULL,
                                        max_year = NULL){
  # Generate species accumulation curves
  #
  # Arguments:
  #   x, numeric, vector of years
  #   ladderise, logical, whether or not to return a ladderised version
  #
  # Returns:
  #   df, data.frame
  
  if(!is.numeric(x)){
    stop("x must be a numeric")
  }
  x <- x[!is.na(x)] # remove any NAs
  
  if(is.null(min_year)){
    min_year <- min(x)
  } else {
    min_year <- min_year
  }
  if(is.null(max_year)){
    max_year <- max(x)
  } else {
    max_year <- max_year
  }
  
  years <- (min_year-1):max_year
  nspp_described <- vector()
  for(i in 1:length(years)){
    nspp_described[i] <- length(x[x <= years[i]])
  }
  
  if(ladderise){
    nspp_step <- rep(nspp_described, each = 2)[-length(nspp_described)*2]
    t_step <- rep(years, each = 2)[-1]
    return(data.frame("S" = nspp_step, "t" = t_step))
    
  } else {
    return(data.frame("S" = nspp_described, "t" = years))
  }
}