####################################################################################################
### Function to read in intersections and do some basic processing:
### - merges with group info from Global_totals.csv
### - removes australia and australia-only species
### - removes unncessary columns
### - merges Andamans to IndoChina, and Palawan to Philippines
###
### output: a list with the mean ('mean'), and for range maps also lower ('lower') and
###         upper ('upper') presence-absence in each bioregion/marine ecoregion
####################################################################################################

processIntersection <- function(taxon, taxonInfo, dataDir) {
  require(readr)
  require(dplyr)
  
  ### Extract taxon info from global info
  taxType <- taxonInfo$Type[taxonInfo$Taxon == taxon]
  taxName <- taxonInfo$Name[taxonInfo$Taxon == taxon]
  
  ### For checklists we just need the single csv file
  if(taxType %in% c("Checklist", "Occurrence points")) {
    res <- read_csv(file.path(data.dir, paste0("Intersections_bioregions_", taxName, ".csv")),
                    progress = FALSE, show_col_types = FALSE) %>%
      select(-any_of(c("AsiaNative", "OutsideAsia")))
  }
  
  ### For range maps we need have mean, upper and lower estimates
  if(taxType == "Range maps") {
    res <- read_csv(file.path(data.dir, paste0("Intersections_bioregions_", taxName, "_5.0_25.csv")),
                       progress = FALSE, show_col_types = FALSE) %>%
      select(-any_of(c("AsiaNative", "OutsideAsia")))
    resLow <- read_csv(file.path(data.dir, paste0("Intersections_bioregions_", taxName, "_7.5_35.csv")),
                       progress = FALSE, show_col_types = FALSE) %>%
      select(-any_of(c("AsiaNative", "OutsideAsia")))
    resUpp <- read_csv(file.path(data.dir, paste0("Intersections_bioregions_", taxName, "_2.5_15.csv")),
                       progress = FALSE, show_col_types = FALSE) %>%
      select(-any_of(c("AsiaNative", "OutsideAsia")))
  }
  
  ### Remove australia and any species only found there
  res <- res %>%
    select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
    filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic")))) > 0)
  
  if(taxType == "Range maps") {
    resLow <- resLow %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic")))) > 0)
    
    resUpp <- resUpp %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic")))) > 0)
  }
  
  ### Merge Palawan into Philippines
  if(any(names(res) == "Palawan")) {
    res <- res %>%
      mutate(Philippines = case_when(Philippines == 0 & Palawan == 0 ~ 0,
                                     Philippines == 1 | Palawan == 1 ~ 1)) %>%
      select(-Palawan)
    
    if(taxType == "Range maps") {
      resLow <- resLow %>%
        mutate(Philippines = case_when(Philippines == 0 & Palawan == 0 ~ 0,
                                       Philippines == 1 | Palawan == 1 ~ 1)) %>%
        select(-Palawan)
      
      resUpp <- resUpp %>%
        mutate(Philippines = case_when(Philippines == 0 & Palawan == 0 ~ 0,
                                       Philippines == 1 | Palawan == 1 ~ 1)) %>%
        select(-Palawan)
    }
  }
  
  ### Merge Andamans into IndoChina
  if(any(names(res) == "Andamans")) {
    res <- res %>%
      mutate(IndoChina = case_when(IndoChina == 0 & Andamans == 0 ~ 0,
                                   IndoChina == 1 | Andamans == 1 ~ 1)) %>%
      select(-Andamans)
    
    if(taxType == "Range maps") {
      resLow <- resLow %>%
        mutate(IndoChina = case_when(IndoChina == 0 & Andamans == 0 ~ 0,
                                     IndoChina == 1 | Andamans == 1 ~ 1)) %>%
        select(-Andamans)
      
      resUpp <- resUpp %>%
        mutate(IndoChina = case_when(IndoChina == 0 & Andamans == 0 ~ 0,
                                     IndoChina == 1 | Andamans == 1 ~ 1)) %>%
        select(-Andamans)
    }
  }
  
  ### Get info on groups etc from global_totals
  res <- res %>%
    mutate(Taxon = taxon) %>%
    left_join(taxonInfo[, c("Group", "Taxon_Higher", "Taxon", "Global_total")], by = "Taxon") %>%
    select(-any_of(c("Range_area"))) %>%
    relocate(Taxon, Group, Taxon_Higher, Global_total, Species, Genus, Family, Year, AsiaEndemic)

  if(taxType == "Range maps") {
    resLow <- resLow %>%
      mutate(Taxon = taxon) %>%
      left_join(taxonInfo[, c("Group", "Taxon_Higher", "Taxon", "Global_total")], by = "Taxon") %>%
      select(-any_of(c("Range_area"))) %>%
      relocate(Taxon, Group, Taxon_Higher, Global_total, Species, Genus, Family, Year, AsiaEndemic)
    
    resUpp <- resUpp %>%
      mutate(Taxon = taxon) %>%
      left_join(taxonInfo[, c("Group", "Taxon_Higher", "Taxon", "Global_total")], by = "Taxon") %>%
      select(-any_of(c("Range_area"))) %>%
      relocate(Taxon, Group, Taxon_Higher, Global_total, Species, Genus, Family, Year, AsiaEndemic)
  }

  if(taxType %in% c("Checklist", "Occurrence points")) {
    return(list("mean" = res))
  }
  if(taxType == "Range maps") {
    return(list("mean"  = res,
                "lower" = resLow,
                "upper" = resUpp))
  }
}

# processIntersection(taxon = taxon, taxonInfo = taxonInfo, dataDir = data.dir)
