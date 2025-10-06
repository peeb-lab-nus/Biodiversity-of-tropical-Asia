####################################################################################################
###
### Takes a string that includes distributional information, and extracts any words that match up
### with a vector of names for countries, provinces etc (including where provinces are inside
### parantheses) and discards all other information. Only tested for trichoptera data so far.
###
### Returns vector of comma-separated names, maintaining parantheses where appropriate
###
### Charlie Marsh
### charliem2003@github
### 07/2024
###
####################################################################################################

extractRegionInfo <- function(distribution,  # character string with distribution information
                              regionNames    # vector of country/province/region names to check against
) {
  require(stringr)
  
  ### make everything lower case and make sure we have versions with & without special characters
  regionNames <- tolower(regionNames)
  regionNames <- unique(c(regionNames, iconv(regionNames, from = 'UTF-8', to = 'ASCII//TRANSLIT')))
  x <- tolower(distribution)
  
  ### find anything that matches any name in the list of country/province/region names
  matches <- regionNames[sapply(regionNames, function(i) grepl(paste0("\\<", i, "\\>"), x, ignore.case = TRUE))]
  matches <- matches[!matches %in% c("na", "north", "south", "east", "west",
                                     "northern", "southern", "eastern", "western")]
  
  if(length(matches) > 0) {
    ### make sure there are commas after names
    for(m in matches) {
      if(sum(grepl(paste0("\\<", m, "\\>"), matches)) == 1) {  # if match is a subset of larger match then skip
        allLocs <- as.data.frame(str_locate_all(x, m))
        for(i in 1:nrow(allLocs)) {
          loc <- allLocs[i, ]
          if((str_sub(x, loc$end + 1, loc$end + 1) == " ")) {#& (str_sub(x, loc$end + 2, loc$end + 2) != "(")) {
            str_sub(x, loc$end + 1, loc$end + 1) <- ", "
          }
        }
      }
    }
    
    ### get locations of all matches
    allLocs <- data.frame()
    for(m in matches) {
      locs <- as.data.frame(str_locate_all(x, m))
      if(nrow(locs) == 0) { next }
      
      ### include before and after punctuation if necessary
      for(i in 1:nrow(locs)) {
        loc <- locs[i, ]
        if(str_sub(x, loc$start - 1, loc$start - 1) %in% c("(", " ")) { loc$start <- loc$start - 1 }
        if(str_sub(x, loc$end + 1, loc$end + 1) %in% c(",", ".", " ")) { loc$end <- loc$end + 1 }
        if(str_sub(x, loc$end + 1, loc$end + 1) %in% c(")")) { loc$end <- loc$end + 2 }
        # print(c(paste(strsplit(x, "")[[1]][sort(unique(loc$start:loc$end))], collapse = ""), loc$start, loc$end))
        allLocs <- rbind(allLocs, loc)
      }
    }
    
    ### Remove overlaps/duplicates and join together
    allLocs <- sort(unique(unlist(apply(allLocs, 1, function(m) m[1]:m[2]))))
    if(str_sub(x, allLocs[1], allLocs[1]) == " ") { allLocs <- allLocs[-1] }
    x <- paste(strsplit(x, "")[[1]][allLocs], collapse = "")
    
    ### Some tidying up of punctuation
    if(str_sub(x[length(x)], nchar(x[length(x)]), nchar(x[length(x)])) == " ") {
      x[length(x)] <- substr(x[length(x)], 1, nchar(x[length(x)]) - 1)
    }
    x <- gsub(",{2,}", ",", x, perl = TRUE)
    x <- gsub(") ", "), ", x, fixed = TRUE)
    x <- gsub(".", ",", x, fixed = TRUE)
    x <- gsub("  ", " ", x, fixed = TRUE)
    x <- gsub(",(", " (", x, fixed = TRUE)
    x <- gsub("NA", "", x)
    if(str_sub(x[length(x)], nchar(x[length(x)]), nchar(x[length(x)])) == ",") {
      x[length(x)] <- substr(x[length(x)], 1, nchar(x[length(x)]) - 1)
    }
    x <- str_to_title(x)
    # x
    return(as.vector(x))
    
    # ### Old method
    # if(length(matches) > 0) {
    #   ### only take words matching those names and rejoin
    #   x <- unlist(strsplit(x, " "))
    #   x <- x[gsub("[[:punct:]]", "", x, perl = TRUE) %in% matches]
    #   x <- paste(x, collapse = " ")
    #   
    #   ### Some tidying up
    #   x <- gsub(",,", ",", x, perl = TRUE)
    #   x <- gsub(") ", "), ", x, fixed = TRUE)
    #   x <- gsub(":", ",", x, fixed = TRUE)
    #   x <- gsub("na, ", "", x)
    #   
    #   ### If necessary add in commas to separate names, but not if there is already a match
    #   x <- unlist(strsplit(x, ", "))
    #   for(i in 1:length(x)) {
    #     matches <- x[i] %in% regionNames
    #     if(matches == FALSE) { 
    #       j <- unlist(strsplit(x[i], " "))
    #       if(length(j) > 1) {
    #         for(k in 1:(length(j) - 1)) {
    #           if(grepl("(", j[k + 1], fixed = TRUE) == FALSE &
    #              grepl("(", j[k], fixed = TRUE) == FALSE &
    #              grepl(")", j[k + 1], fixed = TRUE) == FALSE) { j[k] <- paste0(j[k], ",") }
    #         }
    #       }
    #       x[i] <- paste(j, collapse = " ")
    #     }
    #   }
    #   if(substr(x[length(x)], nchar(x[length(x)]), nchar(x[length(x)])) == ",") {
    #     x[length(x)] <- substr(x[length(x)], 1, nchar(x[length(x)]) - 1)
    #   }
    #   x <- paste(x, collapse = ", ")
    #   x <- tools::toTitleCase(x)
    #   return(x)
  } else {
    return(NA)
  }
}



### Function for extracting regional data from inside parantheses and append to country
### e.g. Indonesia (Sumatra, Java) will become Indonesia Sumatra, Indonesia Java
# 
# extractRegionInfo <- function(distribution) {
#   splits <- strsplit(distribution, ", ")[[1]]
#   splits <- as.vector(sapply(splits, function(x) trimws(x, "both")))
#   if(sum(!is.na(stringr::str_match(splits, "[(]"))) == 
#      sum(!is.na(stringr::str_match(splits, "[)]")))) {
#     brackets <- data.frame(open  = which(!is.na(stringr::str_match(splits, "[(]"))),
#                            close = which(!is.na(stringr::str_match(splits, "[)]"))))
#     if(nrow(brackets) > 0) {
#       for(bracket in 1:nrow(brackets)) {
#         split <- splits[brackets$open[bracket]:brackets$close[bracket]]
#         if(length(split) == 1) {
#           ### if brackets only contain one region then simply remove brackets
#           split <- gsub("[(]|[)]", "", split)
#           splits[brackets$open[bracket]:brackets$close[bracket]] <- split
#         }
#         if(length(split) > 1) {
#           ### if brackets contain more than one region then paste country name to start of each
#           country <- strsplit(split[1], "[ (]")[[1]][1]
#           split1 <- strsplit(split[1], "[ (]")[[1]]
#           split[1] <- paste(split1[3:length(split1)], collapse = " ")
#           split <- gsub("[(]|[)]", "", split)
#           split <- as.vector(sapply(split, function(x) trimws(x, "both")))
#           split <- paste(country, split, sep = " ")
#           splits[brackets$open[bracket]:brackets$close[bracket]] <- split
#         }
#       }
#     }
#   }
#   final <- paste(splits, collapse = ", ")
#   
#   ### some cleaning up
#   final <- gsub('["]', "", final)
#   final <- gsub('[(]', "", final)
#   final <- gsub('[)]', "", final)
#   return(final)
# }
