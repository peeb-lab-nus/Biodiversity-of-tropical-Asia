####################################################################################################
### Set of functions for splitting up range maps into smaller chunks for quicker intersections
### Charlie Marsh (written when part of Jetz lab, Yale)
### charliem2003@github
### 05/2024
###
####################################################################################################

#==================================================================================================#
#--- Function to calculate the area of the range map bounding box ---------------------------------#
#==================================================================================================#

shpBBoxArea <- function(speciesShp) {
  exts <- st_make_valid(st_as_sfc(st_bbox(speciesShp)))
  shpArea <- drop_units(set_units(st_area(exts), km^2))
  return(shpArea)
}

#==================================================================================================#
#--- Function to split the species shapefile in to quarters ---------------------------------------#
#--- Returns a list of four shapefiles                      ---------------------------------------#
#==================================================================================================#

qrtrShp <- function(speciesShp) {
  ### define the quarters
  exts <- st_bbox(speciesShp)
  qrtrs <- list(LTop = c(xmin = as.numeric(exts$xmin), ymin = as.numeric((exts$ymax + exts$ymin) / 2),
                         xmax = as.numeric((exts$xmax + exts$xmin) / 2), ymax = as.numeric(exts$ymax)),
                RTop = c(xmin = as.numeric((exts$xmax + exts$xmin) / 2), ymin = as.numeric((exts$ymax + exts$ymin) / 2),
                         xmax = as.numeric(exts$xmax), ymax = as.numeric(exts$ymax)),
                LBot = c(xmin = as.numeric(exts$xmin), ymin = as.numeric(exts$ymin),
                         xmax = as.numeric((exts$xmax + exts$xmin) / 2), ymax = as.numeric((exts$ymax + exts$ymin) / 2)),
                RBot = c(xmin = as.numeric((exts$xmax + exts$xmin) / 2), ymin = as.numeric(exts$ymin),
                         xmax = as.numeric(exts$xmax), ymax = as.numeric((exts$ymax + exts$ymin) / 2)))
  
  ### split the shapefile
  speciesShp_qrtrs <- lapply(qrtrs, function(x) {
    temp <- speciesShp %>%
      st_crop(st_bbox(x)) %>%
      st_make_valid() 
    return(temp)
  })
  
  ### remove any quarters with no areas
  # qrtrsAreas <- lapply(speciesShp_qrtrs, function(x) shpBBoxArea(x))
  # qrtrsAreas <- lapply(speciesShp_qrtrs, function(x) drop_units(set_units(st_area(x), km^2)))
  # qrtrsAreas <- unlist(as.numeric(lapply(qrtrsAreas, sum)))
  # speciesShp_qrtrs <- speciesShp_qrtrs[which(qrtrsAreas > 0)]
  speciesShp_qrtrs <- speciesShp_qrtrs[which(unlist(lapply(speciesShp_qrtrs, nrow)) > 0)]
  
  return(speciesShp_qrtrs)
}

#==================================================================================================#
#--- Function to run multiple splits recursively (up to 3 iterations) until no one shapefile ------#
#--- exceeds the threshold                                                                   ------#
#==================================================================================================#

splitShp <- function(speciesShp,  areaThresh, verbose = TRUE) {
  ### check size of initial shapefile
  spArea <- shpBBoxArea(speciesShp)
  if(spArea <= areaThresh) {
    if(verbose == TRUE) {
      print(paste("Initial shapefile below threshold of ", areaThresh, "... no splitting required"))
    }
    speciesShp <- list(speciesShp)
  }
  if(spArea > areaThresh) {
    if(verbose == TRUE) {
      print(paste("Initial Shapefile over threshold of ", areaThresh, "... splitting in to quarters"))
    }
    speciesShp <- qrtrShp(speciesShp)
    
    ### check whether we still need splitting
    qrtrsAreas <- unlist(lapply(speciesShp, shpBBoxArea))
    names(qrtrsAreas) <- gsub("[.max]", "", names(qrtrsAreas))
    needSplitting <- names(qrtrsAreas)[which(qrtrsAreas > areaThresh)]
    if(length(needSplitting) > 0) {
      for(qrtr in needSplitting) {
        if(verbose == TRUE) { print(paste("Splitting", qrtr, "quarter further")) }
        
        speciesShp_qrtr <- speciesShp[[qrtr]]
        speciesShp_qrtr <- qrtrShp(speciesShp_qrtr)
        names(speciesShp_qrtr) <- paste(qrtr, names(speciesShp_qrtr), sep = "_")
        
        ### check none of these exceed the threshold
        qrtrsAreas <- unlist(lapply(speciesShp_qrtr, shpBBoxArea))
        names(qrtrsAreas) <- gsub("[.max]", "", names(qrtrsAreas))
        needSplitting_qrtr <- names(qrtrsAreas)[which(qrtrsAreas > areaThresh)]
        if(length(needSplitting_qrtr) > 0) {
          for(qrtr_qrtr in needSplitting_qrtr) {
            if(verbose == TRUE) { print(paste("Splitting", qrtr_qrtr, "quarter even further")) }
            
            speciesShp_qrtr_qrtr <- speciesShp_qrtr[[qrtr_qrtr]]
            speciesShp_qrtr_qrtr <- qrtrShp(speciesShp_qrtr_qrtr)
            names(speciesShp_qrtr_qrtr) <- paste(qrtr_qrtr, names(speciesShp_qrtr_qrtr), sep = "_")
            
            ### replace the quarter that needed splitting with the new quarter
            speciesShp_qrtr <- speciesShp_qrtr[names(speciesShp_qrtr) != qrtr_qrtr]
            speciesShp_qrtr <- c(speciesShp_qrtr, speciesShp_qrtr_qrtr)
          }
        }
        
        ### replace the quarter that needed splitting with the new quarter
        speciesShp <- speciesShp[names(speciesShp) != qrtr]
        speciesShp <- c(speciesShp, speciesShp_qrtr)
      }
    }
  }
  if(verbose == TRUE) { print("Splitting shapefile complete") }
  return(speciesShp)
}
