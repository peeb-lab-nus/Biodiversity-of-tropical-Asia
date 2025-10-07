####################################################################################################
### Cleaning GBIF data for freshwater crabs
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### GBIF data used: occurrence.parquet - official data dump 01/10/2024
### Flags used:     cleaningSteps_BTAS_10_2024.csv
###
### CITATION: GBIF.org (01 October 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.av6eex
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
occDir  <- file.path("GBIF", "data", "dump", "directory")             # dir that contains the GBIF data dump

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

### Make dataset connection with parquet file of full GBIF database (3,043,244,961 rows)
gbif_par <- open_dataset(file.path(occDir, "gbif_2024_10_01", "occurrence.parquet"))

### Filter just the freshwater crabs and remove points outside region (31,695 rows globally)
gbif <- gbif_par %>%
  filter(family %in% c("Deckeniidae",        #    597 records
                       "Gecarcinucidae",     #  4,149 records
                       "Potamidae",          # 12,490 records
                       "Potamonautidae",     #  5,888 records
                       "Pseudothelphusidae", #  5,027 records
                       "Trichodactylidae"    #  3,544 records
  )) %>%
  select(gbifid, species, genus, family, order, class, phylum, kingdom,
         taxonrank, taxonkey, scientificname, occurrencestatus, basisofrecord, typestatus,
         decimallongitude, decimallatitude, coordinateuncertaintyinmeters, coordinateprecision,
         countrycode, locality, stateprovince, elevation, elevationaccuracy,
         institutioncode, collectioncode, issue, eventdate, year, month, day) %>%
  as_tibble()
table(gbif$family)

### Clean dates
gbif <- cleanDates(gbif = gbif)

### Unlist issues and add in new columns for Flags and Step numbers
gbif <- gbif %>%
  mutate(issue = sapply(issue, function(x) paste(x, collapse = ";")[[1]])) %>%
  mutate(Flag = NA) %>%
  mutate(Flag_step = NA)
# # unique(unlist(strsplit(gbif$issue, ";")))

### Column info for common arguments
colID  <- "gbifid"             # Column name containing unique record ID
colLon <- "decimallongitude"   # Column name containing longitude coords
colLat <- "decimallatitude"    # Column name containing latitude coords

### Read in data frame with cleaning steps to be applied
filterSteps <- read.csv(file.path(projDir, "Occurrences", "cleaning", "cleaningSteps_BTAS_10_2024.csv")) %>%
  filter(!is.na(Order)) %>%
  arrange(Order)

#==================================================================================================#
#------------------------------------------ Run cleaning ------------------------------------------#
#==================================================================================================#

cleaned <- runCleaningSteps(gbif             = gbif,        # prepared gbif data frame, with cols Flag and Flag_step added
                            filterSteps      = filterSteps, # data frame with details of flags to apply
                            id_col           = colID,       # Column name containing unique record ID
                            lon_colname      = colLon,      # Column name containing longitude coords
                            lat_colname      = colLat,      # Column name containing latitude coords
                            pathFlaggedData  = NULL,        # path and file name of full gbif data with column indicated if flag applied. If NULL saves nothing
                            pathFilteredData = file.path(resDir, "gbif_filtered_freshwater_crabs.csv"), # path and file name of gbif data with flagged data filtered out.
                            plot = FALSE)

#==================================================================================================#
#--------------------------------------------- Tidy up --------------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
.rs.restartR()

