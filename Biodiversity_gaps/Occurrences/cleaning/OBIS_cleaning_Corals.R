####################################################################################################
### Cleaning OBIS data for Sharks and rays
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### obis data used: occurrence.parquet - official data dump 23/07/2024
### Flags used:     cleaningSteps_OBIS_12_2024.csv
###
### CITATION: OBIS (2024) Ocean Biodiversity Information System. Intergovernmental Oceanographic 
###   Commission of UNESCO. https://obis.org. Downloaded on 31/10/2024
###
### saves csv with filtered data
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(readr)
library(dplyr)
library(stringi)
library(lubridate)
library(arrow)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir <- "Biodiversity_gaps"                                        # project dir
occDir  <- file.path("OBIS", "data", "dump", "directory")             # dir that contains the GBIF data dump

### You shouldn't need to adjust these folders
funDir  <- file.path(projDir, "Analysis_functions", "Point_cleaning") # dir that contains the function scripts
resDir  <- file.path(projDir, "Occurrences", "cleaning", "cleaned")   # dir to save results to

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Functions for data cleaning
source(file.path(funDir, "cleaningFuns.R"))
source(file.path(funDir, "runCleaningSteps.R"))
source(file.path(funDir, "cleanDates.R"))

### Make dataset connection with parquet file of full obis database (128,621,187 rows)
obis_par <- open_dataset(file.path(occDir, "obis_20240723.parquet"))

### Filter just the corals (979,709 rows)
obis <- obis_par %>%
  filter(family %in% c("Acroporidae",      # 202,649 records
                       "Agariciidae",      #  74,203 records
                       "Astrangiidae",     #   1,861 records
                       "Astrocoeniidae",   #  10,007 records
                       "Cladocoridae",     #   2,724 records
                       "Coscinaraeidae",   #     745 records
                       "Diploastraeidae",  #   1,542 records
                       "Euphylliidae",     #   8,423 records
                       "Faviidae",         #  64,886 records
                       "Fungiidae",        #  21,857 records
                       "Helioporidae",     #     859 records
                       "Leptastreidae",    #   8,932 records
                       "Lobophylliidae",   #  14,170 records
                       "Meandrinidae",     #  13,661 records
                       "Merulinidae",      # 156,821 records
                       "Milleporidae",     #  23,872 records
                       "Montastraeidae",   #  29,370 records
                       "Oculinidae",       #   8,680 records
                       "Oulastreidae",     #     243 records
                       "Plerogyridae",     #   1,314 records
                       "Plesiastreidae",   #   1,069 records
                       "Pocilloporidae",   # 100,365 records
                       "Poritidae",        # 156,046 records
                       "Psammocoridae",    #  10,433 records
                       "Rhizangiidae",     #  39,035 records
                       "Tubiporidae",      #   2,550 records
                       "Turbinoliidae"     #   4,228 records
                       ) |
          genus %in% c("Balanophyllia",    #   5,086 records
                       "Duncanopsammia",   #     973 records
                       "Heterocyathus",    #     984 records
                       "Heteropsammia",    #     829 records
                       "Pachyseris",       #   1,844 records
                       "Solenastrea",      #   1,149 records
                       "Turbinaria"        #   8,299 records
  )) %>%
  select(id, species, genus, family, order, class, phylum, kingdom,
         taxonRank, scientificName, occurrenceStatus, basisOfRecord, typeStatus,
         decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, coordinatePrecision,
         countryCode, locality, stateProvince, depth,
         institutionCode, collectionCode, flags, eventDate, year, month, day) %>%
  as_tibble()
table(obis$order)
table(obis$family)
table(obis$genus)

### Clean dates
obis <- cleanDates(gbif = obis)

### Unlist flags and add in new columns for Flags and Step numbers
obis <- obis %>%
  rename(issue = flags) %>%   # to match to gbif name
  mutate(Flag = NA) %>%
  mutate(Flag_step = NA)

### Column info for common arguments
colID  <- "id"                 # Column name containing unique record ID
colLon <- "decimalLongitude"   # Column name containing longitude coords
colLat <- "decimalLatitude"    # Column name containing latitude coords

### Read in data frame with cleaning steps to be applied
filterSteps <- read.csv(file.path(projDir, "Occurrences", "cleaning", "cleaningSteps_OBIS_12_2024.csv")) %>%
  filter(!is.na(Order)) %>%
  arrange(Order)

#==================================================================================================#
#------------------------------------------ Run cleaning ------------------------------------------#
#==================================================================================================#

cleaned <- runCleaningSteps(gbif             = obis,        # prepared obis data frame, with cols Flag and Flag_step added
                            filterSteps      = filterSteps, # data frame with details of flags to apply
                            id_col           = colID,       # Column name containing unique record ID
                            lon_colname      = colLon,      # Column name containing longitude coords
                            lat_colname      = colLat,      # Column name containing latitude coords
                            pathFlaggedData  = NULL,        # path and file name of full obis data with column indicated if flag applied. If NULL saves nothing
                            pathFilteredData = file.path(resDir, "obis_filtered_corals.csv"), # path and file name of obis data with flagged data filtered out.
                            plot = FALSE)

#==================================================================================================#
#--------------------------------------------- Tidy up --------------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
.rs.restartR()
