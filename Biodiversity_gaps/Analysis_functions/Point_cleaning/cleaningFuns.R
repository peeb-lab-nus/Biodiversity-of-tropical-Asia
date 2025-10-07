####################################################################################################
### Individual cleaning functions for flagging GBIF data. Each function is called by
### runCleaningSteps, which in turn reads in the filterSteps spreadsheet to determine which flags to
### apply and in which order to run them in.
### 
### Returns:
###   The original GBIF data frame with columns for whether records were flagged ('Flag'), and if so
###    at which stage ('Flag_step')
###
### Charlie Marsh
### charliem2003@github
### 10/2024
####################################################################################################

### To retrieve all unique entries for entire GBIF dump for a given column:
# allFiles <- list.files("/mnt/Documents/occurrence.parquet/", full.names = TRUE)
# all <- vector()
# pb <- txtProgressBar(1, length(allFiles), style = 3)
# for(i in 1:length(allFiles)) {
#   setTxtProgressBar(pb, i)
#   gbifPart <- read_parquet(allFiles[i], as_data_frame = TRUE, col_select = "taxonrank")
#   cats <- arrow_table(gbifPart) %>% unique() %>% collect()
#   all <- unique(c(all, as.vector(cats)[[1]]))
# }
# all

####################################################################################################
### Incomplete or invalid latitude and longitude coordinates
incompleteCoords <- function(gbif,
                             filterSteps,
                             id_col           = "gbifID",           # Column name containing unique record ID
                             lon_colname      = "decimalLongitude", # Column name containing longitude coords
                             lat_colname      = "decimalLatitude",   # Column name containing latitude coords
                             ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "No_coordinates"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with no/incomplete coordinates
  flagged <- gbifClean %>%
    filter(
      (is.na(get(lon_colname))) |
        (is.na(get(lat_colname))) |
        (grepl(issues, issue))
    ) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- incompleteCoords(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Record type is fossil or a living specimen (e.g. in zoo or botanic garden)
recordType <- function(gbif,
                       filterSteps,
                       id_col = "gbifID",           # Column name containing unique record ID
                       ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Record_type"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records labelled as fossils, or living specimens (in zoos, botanic gardens etc)
  flagged <- gbifClean %>%
    filter(if_all(matches("basisOfRecord", ignore.case = TRUE)) %in% c(
      # "MaterialEntity", "MaterialSample",
      # "Event",
      # "Taxon",
      # "OBSERVATION",
      # "PRESERVED_SPECIMEN", "PreservedSpecimen",
      # "MATERIAL_SAMPLE",
      # "HUMAN_OBSERVATION", "humanobservation", "Human observation", "HumanObservation",
      # "NomenclaturalChecklist",
      # "MACHINE_OBSERVATION", "MachineObservation",
      # "OCCURRENCE", "Occurrence",
      # "MATERIAL_CITATION",
      "FOSSIL_SPECIMEN",
      "LIVING_SPECIMEN", "LivingSpecimen"
    )) %>%
    select(all_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- recordType(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Occurrence status is absent
occurrenceStatus <- function(gbif,
                             filterSteps,
                             id_col = "gbifID",           # Column name containing unique record ID
                             ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Occurrence_status"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records recorded as not-present
  flagged <- gbifClean %>%
    filter(if_all(matches("occurrenceStatus", ignore.case = TRUE)) %in% c(
      # "PRESENT", "present", "Presente", "Q", "Present", "P", "presence"
      # "Occasional", "Frequent", "Abundant", "Common", "Rare",
      # "Historic",
      # "1", "11", "72", "250", "380", "408", "476", "510", "630", "730",
      "ABSENT",
      "0"
    )) %>%
    filter(if_any(matches("absence"), \(absence) absence == TRUE)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- occurrenceStatus(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Record not resolved to species-level
taxonomicRank <- function(gbif,
                          filterSteps,
                          id_col = "gbifID",           # Column name containing unique record ID#
                          ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "No_species_name"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with name not at species-level (or lower)
  flagged <- gbifClean %>%
    filter((!if_all(matches("taxonRank", ignore.case = TRUE)) %in% c(
               "UNRANKED",
               "VARIETY", "Variety",
               "FORM", "Forma",
               "SPECIES", "Species",
               "SUBSPECIES", "Subspecies",
               NA
               # "KINGDOM",
               # "PHYLUM",
               # "CLASS", "Class",
               # "ORDER", "Order",
               # "FAMILY", "Family", "Subfamily",
               # "GENUS", "Genus", "Subgenus"
             )) |
             (is.na(species)) |
             (grepl(issues, issue))
    ) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- taxonomicRank(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Record not resolved to species-level or name not accepted
taxonomicStatus <- function(gbif,
                            filterSteps,
                            id_col = "gbifID",           # Column name containing unique record ID#
                            ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Name_status"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with name not at species-level (or lower)
  flagged <- gbifClean %>%
    filter(
      (if_all(matches("taxonomicStatus", ignore.case = TRUE)) %in% c(
        # "ACCEPTED", "accepted", "Aceptado", "Válido",
        # "DOUBTFUL", "Uncertain",
        # NA,
        "Inválido", "unaccepted", "alternate representation",
        "SYNONYM", "Sinónimos"
      )) |
        (grepl(issues, issue))
    ) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- taxonomicStatus(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Records with 0,0 coordinates
zeroCoordinates <- function(gbif,
                            filterSteps,
                            id_col       = "gbifID",           # Column name containing unique record ID
                            lon_colname  = "decimalLongitude", # Column name containing longitude coords
                            lat_colname  = "decimalLatitude",  # Column name containing latitude coords
                            ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Zero_coordinates"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with no coordinates
  flagged <- gbifClean %>%
    filter(
      (get(lon_colname) == 0 & get(lat_colname) == 0) |
        grepl(issues, issue)
    ) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- zeroCoordinates(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Records with no date (either eventDate or Year)
noDate <- function(gbif,
                   filterSteps,
                   id_col = "gbifID",           # Column name containing unique record ID
                   ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "No_date"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with no eventDate or year
  flagged <- gbifClean %>%
    filter(is.na(yearClean) & is.na(monthClean)) %>%
    # filter(is.na(eventDate) & is.na(year)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- noDate(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### ***OBIS ONLY*** Record not in specified habitat
filterHabitat <- function(gbif,
                          filterSteps,
                          id_col = "gbifID",           # Column name containing unique record ID#
                          ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Not_in_habitat"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Extract parameters
  params <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    strsplit("; ") %>%
    unlist()
  params <- data.frame(Habitat = t(data.frame(strsplit(params, " == "))[1, ]),
                       Keep    = t(data.frame(strsplit(params, " == "))[2, ]), row.names = NULL)
  colnames(params) <- c("Habitat", "Keep")
  
  ### Only need habitats we want to filter out (Keep == FALSE)
  keep  <- params[params$Keep == TRUE, ]
  throw <- params[params$Keep == FALSE, ]
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records not in specified habitats
  flagged <- gbifClean
  if(nrow(keep) > 0) {
    flagged <- filter(flagged, apply(flagged[, keep$Habitat],  1, function(x) all(x == FALSE) | all(is.na(x))))
  }
  if(nrow(throw) > 0) {
    flagged <- filter(flagged, apply(flagged[, throw$Habitat], 1, function(x) any(x == TRUE)))
  }
  flagged <- flagged %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- filterHabitat(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Records with possibly dodgy coordinates
uncertainCoords <- function(gbif,
                            filterSteps,
                            id_col       = "gbifID",           # Column name containing unique record ID
                            ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Uncertain_coordinates"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with potential issues with coordinates
  flagged <- gbifClean %>%
    filter(grepl(issues, issue)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- uncertainCoords(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Records with possibly dodgy scientific name
uncertainName <- function(gbif,
                          filterSteps,
                          id_col       = "gbifID",           # Column name containing unique record ID
                          ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Uncertain_name"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with potential name issues
  flagged <- gbifClean %>%
    filter(grepl(issues, issue)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- uncertainName(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Record with possibly dodgy date
uncertainDate <- function(gbif,
                          filterSteps,
                          id_col       = "gbifID",           # Column name containing unique record ID
                          ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Uncertain_date"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with potential date issues
  flagged <- gbifClean %>%
    filter(grepl(issues, issue)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- uncertainDate(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Records where coordinate uncertainty (in km) is greater than threshold
coordUncertainty <- function(gbif,
                             filterSteps,
                             id_col       = "gbifID",           # Column name containing unique record ID
                             lon_colname  = "decimalLongitude", # Column name containing longitude coords
                             lat_colname  = "decimalLatitude",  # Column name containing latitude coords
                             ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Coordinate_accuracy"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Extract parameters
  params <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    strsplit("; ")
  threshDist <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("Dist = ", "", .) %>%
    as.numeric() * 1000         # values in gbif are in m - convert from km
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with high coordinate uncertainty
  flagged <- gbifClean %>%
    filter(if_all(matches("coordinateUncertaintyInMeters", ignore.case = TRUE), ~ .x > threshDist)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- coordUncertainty(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Records where year collected is older than threshold
recordAge <- function(gbif,
                      filterSteps,
                      id_col       = "gbifID",           # Column name containing unique record ID
                      lon_colname  = "decimalLongitude", # Column name containing longitude coords
                      lat_colname  = "decimalLatitude",  # Column name containing latitude coords
                      ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "Record_age"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Extract parameters
  params <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    strsplit("; ")
  threshYear <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("Year = ", "", .) %>%
    as.numeric()
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with age older than threshold
  flagged <- gbifClean %>%
    filter(if_all(matches("yearClean", ignore.case = TRUE), ~ .x < threshYear)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- recordAge(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Points at country centroids - uses CoordinateCleaner
ptsCentroids <- function(gbif,
                         filterSteps,
                         id_col           = "gbifID",           # Column name containing unique record ID
                         lon_colname      = "decimalLongitude", # Column name containing longitude coords
                         lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                         # centroids_rad    = 1000,               # Threshold from centroid (in m)
                         # centroids_detail = "country",          # Test for 'country', 'provinces' or 'both'
                         ...
) {
  require(dplyr)
  require(CoordinateCleaner)
  
  ### Extract the filter step information
  stepFlag   <- "Points_at_country_centroids"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Extract parameters
  params <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    strsplit("; ") %>%
    unlist()
  threshDist <- params[grepl("Dist", params)] %>%
    gsub("Dist = ", "", .) %>%
    as.numeric() * 1000         # cc takes distance as m - convert from km
  centroidType <- params[grepl("centroids_detail", params)] %>%
    gsub("centroids_detail = ", "", .)
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Run coordinate cleaner
  gbifClean$passed <- cc_cen(x      = gbifClean,
                             lon    = lon_colname,
                             lat    = lat_colname,
                             geod   = TRUE,
                             buffer = threshDist,
                             test   = centroidType,
                             value  = "flagged",
                             verbose = FALSE)
  
  ### Get IDs of records that didn't pass
  flagged <- unlist(as.vector(gbifClean[gbifClean$passed == FALSE, id_col]))
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- ptsCentroids(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Points at GBIF headquarters - uses CoordinateCleaner
ptsGBIF <- function(gbif,
                    filterSteps,
                    id_col           = "gbifID",           # Column name containing unique record ID
                    lon_colname      = "decimalLongitude", # Column name containing longitude coords
                    lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                    # centroids_rad    = 1000,               # Threshold from centroid (in m)
                    ...
) {
  require(dplyr)
  require(CoordinateCleaner)
  
  ### Extract the filter step information
  stepFlag   <- "Points_at_GBIF_headquarters"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  threshDist <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("Dist = ", "", .) %>%
    as.numeric() * 1000         # cc takes distance as m - convert from km
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Run coordinate cleaner
  gbifClean$passed <- cc_cen(x      = gbifClean,
                             lon    = lon_colname,
                             lat    = lat_colname,
                             geod   = TRUE,
                             buffer = threshDist,
                             value  = "flagged",
                             verbose = FALSE)
  
  ### Get IDs of records that didn't pass
  flagged <- unlist(as.vector(gbifClean[gbifClean$passed == FALSE, id_col]))
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- ptsGBIF(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Points near institutions - uses CoordinateCleaner
ptsInstitutions <- function(gbif,
                            filterSteps,
                            id_col           = "gbifID",           # Column name containing unique record ID
                            lon_colname      = "decimalLongitude", # Column name containing longitude coords
                            lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                            # centroids_rad    = 1000,               # Threshold from centroid (in m)
                            ...
) {
  require(dplyr)
  require(CoordinateCleaner)
  
  ### Extract the filter step information
  stepFlag   <- "Points_near_institutions"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  threshDist <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("Dist = ", "", .) %>%
    as.numeric() * 1000         # cc takes distance as m - convert from km
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Run coordinate cleaner
  gbifClean$passed <- cc_inst(x      = gbifClean,
                              lon    = lon_colname,
                              lat    = lat_colname,
                              geod   = TRUE,
                              buffer = threshDist,
                              value  = "flagged",
                              verbose = FALSE)
  
  ### Get IDs of records that didn't pass
  flagged <- unlist(as.vector(gbifClean[gbifClean$passed == FALSE, id_col]))
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- ptsInstitutions(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Remove records with duplicate date and coordinates (retains first record in duplicate set)
spatioTempDuplicates <- function(gbif,
                                 filterSteps,
                                 id_col           = "gbifID",           # Column name containing unique record ID
                                 lon_colname      = "decimalLongitude", # Column name containing longitude coords
                                 lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                                 ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag <- "Spatiotemporal_duplicates"
  stepNo   <- filterSteps$Order[filterSteps$Flag == stepFlag]
  tempRes  <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("tempRes = ", "", .) %>%
    tolower()
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### year, month and day must be same (if day or month = NA then not considered duplicate)
  if(tempRes == "day") {
    gbifClean <- gbifClean %>%
      filter(!is.na(yearClean) & !is.na(monthClean) & !is.na(dayClean)) %>%
      mutate(dateDup = paste(yearClean, monthClean, dayClean, sep = "/"))
  }
  
  ### year, month must be same (if month = NA then not considered duplicate)
  if(tempRes == "month") {
    gbifClean <- gbifClean %>%
      filter(!is.na(yearClean) & !is.na(monthClean)) %>%
      mutate(dateDup = paste(yearClean, monthClean, sep = "/"))
  }
  
  ### Only year needs to be the same
  if(tempRes == "year") {
    gbifClean <- gbifClean %>%
      filter(!is.na(yearClean)) %>%
      mutate(dateDup = yearClean)
  }
  
  ### Get IDs of records that are spatio-temporal duplicates
  flagged <- gbifClean %>%
    select(matches(c(id_col, "species", lon_colname, lat_colname, "dateDup"), ignore.case = TRUE)) %>%
    filter(duplicated(.[, names(.) != id_col])) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- spatioTempDuplicates(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Remove records with duplicate coordinates (retains first record in duplicate set)
spatialDuplicates <- function(gbif,
                              filterSteps,
                              id_col           = "gbifID",           # Column name containing unique record ID
                              lon_colname      = "decimalLongitude", # Column name containing longitude coords
                              lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                              ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag <- "Spatial_duplicates"
  stepNo   <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records that are spatio-temporal duplicates
  flagged <- gbifClean %>%
    select(matches(c(id_col, "species", lon_colname, lat_colname), ignore.case = TRUE)) %>%
    filter(duplicated(.[, names(.) != id_col])) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- spatialDuplicates(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### Remove points that fall outside region (bounding box of shapefile)
outsideRegion <- function(gbif,
                          filterSteps,
                          id_col       = "gbifID",           # Column name containing unique record ID
                          lon_colname  = "decimalLongitude", # Column name containing longitude coords
                          lat_colname  = "decimalLatitude",  # Column name containing latitude coords
                          ...
) {
  require(dplyr)
  require(sf)
  
  ### Extract the filter step information
  stepFlag <- "Points_outside_region"
  stepNo   <- filterSteps$Order[filterSteps$Flag == stepFlag]
  shpPath  <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    gsub("shpFile = ", "", .)
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Read in shapefile, take bounding box and transform projection if necessary
  region <- st_read(shpPath) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_transform(4326) %>%
    st_make_valid() %>%
    st_as_sf()
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  # ### Convert remaining points to shapefile (presumed crs = 4326)
  # gbifPts <- gbifClean %>%
  #   select(all_of(c(id_col, lon_colname, lat_colname))) %>%
  #   st_as_sf(coords = c(lon_colname, lat_colname),
  #            crs = 4326)
  # 
  # ### Intersections between points and bbox
  # gbifPts <- st_join(gbifPts, region, join = st_within, left = FALSE)
  # 
  # ### Get IDs of records that didn't pass
  # passed  <- unlist(as.vector(st_drop_geometry(gbifPts)[, id_col]))
  # flagged <- unlist(as.vector(gbifClean[, id_col]))
  # flagged <- flagged[!flagged %in% passed]
  
  ### Intersections between points and bbox
  # inside <- st_join(region, gbifPts, join = st_contains, left = FALSE)
  region <- st_bbox(region)
  inside <- gbifClean %>%
    filter(get(lon_colname) >= region$xmin) %>%
    filter(get(lon_colname) <= region$xmax) %>%
    filter(get(lat_colname) >= region$ymin) %>%
    filter(get(lat_colname) <= region$ymax)
  
  ### IDs of points to flag
  passed  <- unlist(as.vector(st_drop_geometry(inside)[, id_col]))
  flagged <- unlist(as.vector(gbifClean[, id_col]))
  flagged <- flagged[!flagged %in% passed]
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- outsideRegion(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)

####################################################################################################
### ***OBIS ONLY*** Filter records that are identified as being on land according to OBIS
obisOnLand <- function(gbif,
                       filterSteps,
                       id_col = "gbifID",           # Column name containing unique record ID#
                       ...
){
  require(dplyr)
  
  ### Extract the filter step information
  stepFlag   <- "OBIS_on_land"
  stepNo     <- filterSteps$Order[filterSteps$Flag == stepFlag]
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Extract parameters
  issues     <- filterSteps$Issue_flags_applied[filterSteps$Flag == stepFlag]
  issues     <- paste(unlist(strsplit(issues, "; ")), collapse = "|")
  print(paste0(stepNo, ": Flagging records by ", stepFlag, "..."))
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Get IDs of records with potential date issues
  flagged <- gbifClean %>%
    filter(grepl(issues, issue)) %>%
    select(any_of(id_col)) %>%
    as.vector(.) %>%
    unlist()
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- obisOnLand(gbif = gbif, filterSteps = filterSteps, id_col = colID)

####################################################################################################
### Remove points that fall outside land areas
maskToLand <- function(gbif,
                       filterSteps,
                       id_col           = "gbifID",           # Column name containing unique record ID
                       lon_colname      = "decimalLongitude", # Column name containing longitude coords
                       lat_colname      = "decimalLatitude",  # Column name containing latitude coords
                       ...
) {
  require(dplyr)
  require(sf)
  
  ### Extract the filter step information
  stepFlag <- "Non-land"
  stepNo   <- filterSteps$Order[filterSteps$Flag == stepFlag]
  params <- filterSteps$Params[filterSteps$Flag == stepFlag] %>%
    strsplit("; ") %>%
    unlist()
  invert <- params[grepl("inverse", params)] %>%
    gsub("inverse = ", "", .)
  shpPath  <- params[grepl("shpFile", params)] %>%
    gsub("shpFile = ", "", .)
  print(paste0(stepNo, ": Flagging records with ", stepFlag, "..."))
  
  ### Read in shapefile
  land <- st_read(shpPath)
  
  ### Filter out points that are already flagged
  gbifClean <- filter(gbif, is.na(Flag))
  
  ### Convert remaining points to shapefile (presumed crs = 4326)
  gbifPts <- gbifClean %>%
    select(all_of(c(id_col, lon_colname, lat_colname))) %>%
    st_as_sf(coords = c(lon_colname, lat_colname),
             crs = 4326)
  
  ### Transform projection if necessary
  gbifPts <- st_transform(gbifPts, st_crs(land))
  
  ### Intersections between points and land
  inside <- st_join(land, gbifPts, join = st_contains, left = FALSE)
  passed  <- unlist(as.vector(st_drop_geometry(inside)[, id_col]))
  
  ### Get IDs of records that didn't pass (or if inverse = TRUE then records on land that did pass)
  if(invert == FALSE) {
    flagged <- unlist(as.vector(gbifClean[, id_col]))
    flagged <- flagged[!flagged %in% passed]
  } else {
    flagged <- passed
  }
  
  ### Assign flags to original data frame
  gbif <- gbif %>%
    mutate(Flag = case_when(get(id_col) %in% flagged ~ stepFlag,
                            .default = Flag)) %>%
    mutate(Flag_step = case_when(get(id_col) %in% flagged ~ stepNo,
                                 .default = Flag_step))
  print(table(gbif$Flag, useNA = "always"))
  
  ### Return full data frame
  return(gbif)
}

# gbif <- maskToLand(gbif = gbif, filterSteps = filterSteps, id_col = colID,
#                          lon_colname = colLon, lat_colname = colLat)