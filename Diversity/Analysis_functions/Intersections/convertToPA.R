####################################################################################################
### Function to take intersections of range maps and bioregions and convert to presence-absence
### Charlie Marsh
### charliem2003@github
### 05/2024
###
### returns data frame with presence-absence for each bioregion
####################################################################################################

convertToPA <- function(rangeIntersects = intersections, # data frame of occupied area for species range maps x bioregions
                        threshAsia      = NULL,          # NULL or prop of species rangemap that needs to occur within trop. asia to be considered
                        threshBioregion = 0.01,          # prop of bioregion area that needs to be occupied to be considered present in bioregion
                        threshPropRange = 0.25,          # if < threshBioregion, prop of range to occur within bioregion to still be considered present
                        threshAsiaEnd   = 0.99,          # prop. total species range that occurs in tropical asia to be considered regional endemic (Australia is removed)
                        rangeAreaCol    = "Range_area",  # col name that contains species total range areas
                        regionNameCol   = "Bioregion",   # column name where unique region names are kept
                        regions         = regions,       # shapefile with bioregions
                        verbose         = TRUE)          # prints info on number of species retained/removed
{
  require(sf)
  require(dplyr)
  require(tidyr)
  require(units)

  ### Vector of bioregions sampled, plus version without australian outgroups  
  regs <- regions %>%
    st_drop_geometry() %>%
    select(all_of(regionNameCol))
  regs <- as.vector(regs[, 1])
  
  ausRegions <- c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf")
  regsNoAus <- regs[!regs %in% ausRegions]
  
  if(verbose == TRUE) {
    print("Cleaning up intersections...")
    print(paste(sum(rowSums(rangeIntersects[, names(rangeIntersects) %in% regsNoAus], na.rm = TRUE) > 0),
                "species with range maps intersecting tropical Asia excluding N. Australia,",
                nrow(rangeIntersects), "including N. Australia"))
  }
  
  ##################################################################################################
  ### Calculate endemicity and if requested filter where proportion range in Asia < AsiaThresh
  rangeIntersects <- rangeIntersects %>%
    as_tibble() %>%
    mutate(totalAsia = rowSums(across(any_of(regs)), na.rm = TRUE)) %>%
    mutate(totalAsia = totalAsia / get(rangeAreaCol)) %>%
    mutate(totalAsiaNoAus = rowSums(across(any_of(regsNoAus)), na.rm = TRUE)) %>%
    mutate(totalAsiaNoAus = totalAsiaNoAus / get(rangeAreaCol))
  
  if(!is.null(threshAsia)) {
    ### So we can keep species that fail asia threshold but do occur in australia outgroup
    rangeIntersects <- rangeIntersects %>%
      mutate(totalAus = rowSums(across(any_of(ausRegions)), na.rm = TRUE))
    
    if(verbose == TRUE) {
      print(paste(sum((rangeIntersects$totalAsiaNoAus < threshAsia) &
                        (rangeIntersects$totalAus == 0)), "species have <", threshAsia,
                  "total range within tropical Asia and also not present in N. Australia. Removing..."))
    }
    
    ### filter out species
    rangeIntersects <- rangeIntersects %>%
      filter((totalAsiaNoAus >= threshAsia) | (totalAus > 0))
    
    if(verbose == TRUE) {
      print(paste(nrow(rangeIntersects), "species remaining"))
    }
    rangeIntersects <- select(rangeIntersects, -totalAus)
  }
  
  rangeIntersects <- rangeIntersects %>%
    mutate(AsiaEndemic = case_when(totalAsiaNoAus >= threshAsiaEnd ~ 1,
                                   totalAsiaNoAus <  threshAsiaEnd ~ 0)) %>%
    select(-totalAsia, -totalAsiaNoAus)
  
  ##################################################################################################
  ### Second, set area to 0 where proportion of bioregion occupied by each species < bioregionThresh
  
  ### areas for each Bioregion
  regionsArea <- regions %>%
    mutate(BioregionArea = drop_units(set_units(st_area(regions), "km^2"))) %>%
    st_drop_geometry() %>%
    select(all_of(regionNameCol), BioregionArea)
  
  ### Prop. bioregion occupied
  propBioregion <- rangeIntersects
  for(bioregion in names(propBioregion)[names(propBioregion) %in% regs]) {
    propBioregion[, bioregion] <- propBioregion[, bioregion] / regionsArea$BioregionArea[regs == bioregion]
  }
  propBioregion <- propBioregion %>%
    mutate(across(any_of(regs), ~ case_when(. == 0 ~ NA,
                                            . >  0 ~ .))) %>%
    mutate(across(any_of(regs), ~ case_when(. >= threshBioregion  ~ 2,
                                            . <  threshBioregion  ~ 1,
                                            is.na(.)              ~ 0))) %>%
    select(Species, any_of(regs))
  
  ### Prop. total range to occur in each bioregion
  propRange <- rangeIntersects %>%
    mutate(across(any_of(regs), ~ . / rangeIntersects[, rangeAreaCol])) %>%
    mutate(across(any_of(regs), ~ case_when(. == 0 ~ NA,
                                            . >  0 ~ .))) %>%
    mutate(across(any_of(regs), ~ case_when(. >= threshPropRange            ~ 2,
                                            (. > 0) & (. < threshPropRange) ~ 1,
                                            is.na(.)                        ~ 0))) %>%
    select(Species, any_of(regs))
  
  both <- tibble(select(propRange, Species),
                 select(propRange, -Species) + select(propBioregion, -Species)) %>%
    mutate(across(any_of(regs), ~ case_when(.x == 4 ~ 1, # >1% of region occupied & >25% of range in region
                                            .x == 3 ~ 1, # either >1% of region occupied or >25% of range in region
                                            .x == 2 ~ 0, # <1% of region occupied & <25% of range in region
                                            .x == 0 ~ NA # not in region 
    ))) %>%
    mutate(ToRemove = rowSums(across(any_of(regs)), na.rm = TRUE)) %>%
    mutate(ToRemove = case_when(ToRemove >= 1 ~ 0,   # at least one bioregion occupied
                                ToRemove == 0 ~ 1))  # no bioregions occupied
  
  if(verbose == TRUE) {
    print(paste0(sum(both[, names(both) %in% regs] == 0, na.rm = TRUE),
                 " cases where prop. bioregion occupied is <", threshBioregion, " and also <", threshPropRange,
                 " of the species range occurs within that bioregion. Removing those intersections..."))
  }
  
  ### remove rangeIntersects which have value of 1 for both propBioregion and propRange
  rangeIntersects <- rangeIntersects  %>%
    select(-any_of(regs)) %>%
    left_join(both, by = "Species") %>%
    mutate(across(any_of(regs), ~ case_when(is.na(.x) ~ 0,
                                            .x == 0   ~ 0,
                                            .x == 1   ~ 1)))
  
  ### Remove any species where all rangeIntersects were < threshBioregion and < threshPropRange
  if(verbose == TRUE) {
    print(paste0(sum(both$ToRemove == 1), " species where all bioregions failed those criteria.",
                 " No occupied bioregion - removing these species..."))
  }
  rangeIntersects <- rangeIntersects %>%
    filter(ToRemove == 0) %>%
    select(-ToRemove)
  
  ### total number of species in tropical asia
  if(verbose == TRUE) {
    print(paste(sum(rowSums(rangeIntersects[, names(rangeIntersects) %in% regsNoAus]) > 0),
                "species for tropical Asia excluding N. Australia,",
                nrow(rangeIntersects), "species including N. Australia"))
    
    ### number species endemic to area
    print(paste(sum(rangeIntersects$AsiaEndemic),
                "species endemic to tropical Asia (excluding N. Australia)"))
  }
  return(rangeIntersects)
}
