####################################################################################################
### Cleaning GBIF data for freshwater fish
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

### Filter all fish - we have to include even marine groups at this stage (68,666,102 rows globally)
gbif <- gbif_par %>%
  filter(order %in% c("Acipenseriformes",     #     20,074 records - sturgeons, yes
                      "Albuliformes",         #     25,970 records - bonefish
                      "Amiiformes",           #     18,547 records - bowfins, yes
                      "Anguilliformes",       #  1,656,064 records - eels, yes
                      "Ateleopodiformes",     #      2,029 records - tadpole fish
                      "Atheriniformes",       #    357,187 records - yes
                      "Aulopiformes",         #    346,095 records - lizardfish
                      "Batrachoidiformes",    #     60,644 records - toad/frogfish
                      "Beloniformes",         #    253,844 records - both
                      "Beryciformes",         #  1,026,617 records - deepwater
                      "Cetomimiformes",       #      2,454 records - whalefish
                      "Characiformes",        #    937,317 records - yes
                      "Clupeiformes",         #  3,701,042 records - herring
                      "Cypriniformes",        # 11,021,049 records - yes
                      "Cyprinodontiformes",   #    722,069 records - yes
                      "Elopiformes",          #    319,119 records - tarpons
                      "Esociformes",          #    704,861 records - pike
                      "Gadiformes",           #  2,885,309 records - cod
                      "Gasterosteiformes",    #    399,167 records - sticklebacks
                      "Gobiesociformes",      #     32,768 records - blennies
                      "Gonorynchiformes",     #      8,138 records - both
                      "Gymnotiformes",        #     55,600 records - knifefish
                      "Lampriformes",         #     10,718 records - ribbonfish
                      "Lepisosteiformes",     #    172,249 records - gars
                      "Lophiiformes",         #    174,665 records - anglerfish
                      "Mugiliformes",         #    541,310 records - mullets
                      "Myctophiformes",       #    360,612 records - deep sea
                      "Notacanthiformes",     #     34,295 records - deep sea
                      "Ophidiiformes",        #    170,529 records - cusk-eels
                      "Osmeriformes",         #    706,751 records - yes
                      "Osteoglossiformes",    #     44,197 records - yes
                      "Perciformes",          # 23,326,011 records - perches, groupers
                      "Percopsiformes",       #     35,599 records - trout perch, yes
                      "Pleuronectiformes",    #  2,763,339 records - flatfish
                      "Polymixiiformes",      #      4,044 records - deep sea
                      "Polypteriformes",      #      4,702 records - birchirs, yes     
                      "Saccopharyngiformes",  #      5,749 records - deep sea
                      "Salmoniformes",        #  5,535,703 records - salmonids
                      "Scorpaeniformes",      #  3,389,650 records - scorpion fish     
                      "Siluriformes",         #  1,424,736 records - catfish, yes     
                      "Stephanoberyciformes", #     26,131 records - deep sea
                      "Stomiiformes",         #    220,635 records - deep sea        
                      "Synbranchiformes",     #     21,129 records - swamp eels
                      "Syngnathiformes",      #    319,436 records - sea horses
                      "Tetraodontiformes",    #  1,219,023 records - pufferfish
                      "Zeiformes"             #     72,557 records - dories
  ) |  class %in% c("Coelacanthi",            #      3,782 records - coelacanths
                    "Dipneusti",              #     12,576 records - lungfish
                    "Elasmobranchii",         #  3,041,647 records - sharks and rays
                    "Holocephali",            #    170,575 records - rabbitfish
                    "Myxini",                 #     18,990 records - hagfish
                    "Petromyzonti"            #    278,798 records - lampreys (NA = 65,139,734)
  )) %>%
  select(gbifid, species, genus, family, order, class, phylum, kingdom,
         taxonrank, taxonkey, scientificname, occurrencestatus, basisofrecord, typestatus,
         decimallongitude, decimallatitude, coordinateuncertaintyinmeters, coordinateprecision,
         countrycode, locality, stateprovince, elevation, elevationaccuracy,
         institutioncode, collectioncode, issue, eventdate, year, month, day) %>%
  as_tibble()
table(gbif$class, useNA = "always")
table(gbif$order[is.na(gbif$class)])
nrow(gbif)

### Too many rows to process. Remove points definitely not in region (4,791,155 records in subset)
gbif <- gbif %>%
  filter(between(decimallongitude,  60, 160)) %>%
  filter(between(decimallatitude,  -25,  40)) %>%
  as_tibble()
nrow(gbif)

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
                            pathFilteredData = file.path(resDir, "gbif_filtered_freshwater_fish.csv"), # path and file name of gbif data with flagged data filtered out.
                            plot = FALSE)

#==================================================================================================#
#--------------------------------------------- Tidy up --------------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
.rs.restartR()
