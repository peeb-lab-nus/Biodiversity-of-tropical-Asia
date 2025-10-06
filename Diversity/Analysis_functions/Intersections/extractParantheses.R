####################################################################################################
###
### Function for extracting region/province data from inside parantheses and append to country name
### e.g. 'Indonesia (Sumatra, Java)' will become 'Indonesia Sumatra, Indonesia Java'
###
### Charlie Marsh
### charliem2003@github
### 05/2024
###
####################################################################################################

extractParantheses <- function(distribution) {
  splits <- strsplit(distribution, ", ")[[1]]
  splits <- as.vector(sapply(splits, function(x) trimws(x, "both")))
  if(sum(!is.na(stringr::str_match(splits, "[(]"))) ==
     sum(!is.na(stringr::str_match(splits, "[)]")))) {
    brackets <- data.frame(open  = which(!is.na(stringr::str_match(splits, "[(]"))),
                           close = which(!is.na(stringr::str_match(splits, "[)]"))))
    if(nrow(brackets) > 0) {
      for(bracket in 1:nrow(brackets)) {
        split <- splits[brackets$open[bracket]:brackets$close[bracket]]
        if(length(split) == 1) {
          ### if brackets only contain one region then simply remove brackets
          split <- gsub("[(]|[)]", "", split)
          splits[brackets$open[bracket]:brackets$close[bracket]] <- split
        }
        if(length(split) > 1) {
          ### if brackets contain more than one region then paste country name to start of each
          country <- strsplit(split[1], "[ (]")[[1]][1]
          split1 <- strsplit(split[1], "[ (]")[[1]]
          split[1] <- paste(split1[3:length(split1)], collapse = " ")
          split <- gsub("[(]|[)]", "", split)
          split <- as.vector(sapply(split, function(x) trimws(x, "both")))
          split <- paste(country, split, sep = " ")
          splits[brackets$open[bracket]:brackets$close[bracket]] <- split
        }
      }
    }
  }
  final <- paste(splits, collapse = ", ")

  ### some cleaning up
  final <- gsub('["]', "", final)
  final <- gsub('[(]', "", final)
  final <- gsub('[)]', "", final)
  return(final)
}
