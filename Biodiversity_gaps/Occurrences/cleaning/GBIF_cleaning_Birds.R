####################################################################################################
### Cleaning GBIF data for Birds
### This version cleans in parallel as there are two many points to clean in one go
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
require(dplyr)
library(lubridate)
library(arrow)
library(parallel)

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

### The bird data is so huge (2 billion points) we will save subsets in feather format
### Generate new folder to save subsets to
if(!dir.exists(file.path(resDir, "bird_subsets_orig"))) {
  dir.create(file.path(resDir, "bird_subsets_orig"))
}

### Filter just the birds and remove points outside region (2,000,742,086 rows globally, 105,173,220 records in region)
gbif_par %>%
  filter(class == "Aves") %>%
  select(gbifid, species, genus, family, order, class, phylum, kingdom,
         taxonrank, taxonkey, scientificname, occurrencestatus, basisofrecord, typestatus,
         decimallongitude, decimallatitude, coordinateuncertaintyinmeters, coordinateprecision,
         countrycode, locality, stateprovince, elevation, elevationaccuracy,
         institutioncode, collectioncode, issue, eventdate, year, month, day) %>%
  filter(between(decimallongitude,  60, 160)) %>%
  filter(between(decimallatitude,  -25,  40)) %>%
  
  ### Save to feather format in 10,000 row chunks
  write_dataset(file.path(resDir, "bird_subsets_orig"), format = "feather", max_rows_per_file = 10000)

#==================================================================================================#
#------------------------------------ Run cleaning in parallel ------------------------------------#
#==================================================================================================#

### Read in data frame with cleaning steps to be applied
filterSteps <- read.csv(file.path(projDir, "Occurrences", "cleaning", "cleaningSteps_BTAS_10_2024.csv")) %>%
  filter(!is.na(Order)) %>%
  arrange(Order)

### Column info for common arguments
colID  <- "gbifid"             # Column name containing unique record ID
colLon <- "decimallongitude"   # Column name containing longitude coords
colLat <- "decimallatitude"    # Column name containing latitude coords

### Number of cores to run in parallel - NOTE: on windows mcapply will only run on 1 thread
ncores <- 18

### List of all the subset feather files
allSubsets <- list.files(file.path(resDir, "bird_subsets_orig"), pattern = "part-")
allSubsets <- gsub(".arrow", "", allSubsets)

### Run parallel cleaning
mclapply(as.list(allSubsets),
         function(i) {
           ### Skip if already done
           outputFile <- file.path(resDir, "bird_subsets_cleaned", paste0(i, ".csv"))
           
           if(!file.exists(outputFile)) {
             ### read in chunk of gbif subset
             gbifSubset <- read_feather(file.path(resDir, "bird_subsets_orig", paste0(i, ".arrow")))
             
             ### Clean dates
             gbifSubset <- cleanDates(gbif = gbifSubset)
             
             ### Unlist issues and add in new columns for Flags and Step numbers
             gbifSubset <- gbifSubset %>%
               mutate(issue = sapply(issue, function(x) paste(x, collapse = ";")[[1]])) %>%
               mutate(Flag = NA) %>%
               mutate(Flag_step = NA)
             
             cleaned <- runCleaningSteps(gbif             = gbifSubset,  # prepared gbif data frame, with cols Flag and Flag_step added
                                         filterSteps      = filterSteps, # data frame with details of flags to apply
                                         id_col           = colID,       # Column name containing unique record ID
                                         lon_colname      = colLon,      # Column name containing longitude coords
                                         lat_colname      = colLat,      # Column name containing latitude coords
                                         pathFlaggedData  = NULL,        # path and file name of full gbif data with column indicated if flag applied. If NULL saves nothing
                                         pathFilteredData = outputFile,  # path and file name of gbif data with flagged data filtered out.
                                         plot = FALSE)
           }
           return(NULL)
         },
         mc.cores = ncores)

#==================================================================================================#
#---------------------------------------- Collate subsets -----------------------------------------#
#==================================================================================================#

### List of subsets
subsets <- list.files(file.path(resDir, "bird_subsets"), pattern = "gbif_filtered_birds", full.names = TRUE)
subsets <- subsets[!grepl("gbif_filtered_birds_subset.csv", subsets)]

### Make sure the column types are correct
colTypes <- read_csv(file.path(resDir, "bird_subsets", "gbif_filtered_birds_subset.csv"), n_max = 100000) %>%
  summarise_all(class) %>%
  as.vector() %>%
  unlist()

### Append together
allSubs <- do.call("bind_rows",
                   lapply(list(1:length(subsets)),
                          function(x) {
                            return(read_csv(subsets[x], col_types = colTypes,
                                            progress = FALSE, show_col_types = FALSE))
                          }
                   ))

### Now all records are together, we have to do one last check for spatio-temporal duplicates
filterStepsDup <- filterSteps %>%
  filter(Flag == "Spatiotemporal_duplicates")

cleaned <- runCleaningSteps(gbif             = allSubs,         # prepared gbif data frame, with cols Flag and Flag_step added
                            filterSteps      = filterStepsDup,  # data frame with details of flags to apply
                            id_col           = colID,           # Column name containing unique record ID
                            lon_colname      = colLon,          # Column name containing longitude coords
                            lat_colname      = colLat,          # Column name containing latitude coords
                            pathFlaggedData  = NULL,            # path and file name of full gbif data with column indicated if flag applied. If NULL saves nothing
                            pathFilteredData = file.path(resDir, "gbif_filtered_birds.csv"),  # path and file name of gbif data with flagged data filtered out.
                            plot = FALSE)

#==================================================================================================#
#--------------------------------------------- Tidy up --------------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
.rs.restartR()
