####################################################################################################
### Cleaning GBIF data for Angiosperms
### This version cleans in parallel
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

### Filter just the angiosperms and remove points outside region (447,723,919 rows globally)
gbif <- gbif_par %>%
  filter(class %in% c("Magnoliopsida",  # 346,893,918 records
                      "Liliopsida"      # 100,830,001 records
  )) %>%
  select(gbifid, species, genus, family, order, class, phylum, kingdom,
         taxonrank, taxonkey, scientificname, occurrencestatus, basisofrecord, typestatus,
         decimallongitude, decimallatitude, coordinateuncertaintyinmeters, coordinateprecision,
         countrycode, locality, stateprovince, elevation, elevationaccuracy,
         institutioncode, collectioncode, issue, eventdate, year, month, day) %>%
  filter(between(decimallongitude,  60, 160)) %>%
  filter(between(decimallatitude,  -25,  40)) %>%
  as_tibble()
table(gbif$class)
nrow(gbif)
# Liliopsida  Magnoliopsida
#  4,337,049     10,870,086   15,207,135 records in regional subset

### Clean dates
gbif <- cleanDates(gbif = gbif)

### Unlist issues and add in new columns for Flags and Step numbers
gbif <- gbif %>%
  mutate(issue = sapply(issue, function(x) paste(x, collapse = ";")[[1]])) %>%
  mutate(Flag = NA) %>%
  mutate(Flag_step = NA)
# # unique(unlist(strsplit(gbif$issue, ";")))

### Generate new folder to save subsets to
if(!dir.exists(file.path(resDir, "angiosperm_subsets"))) {
  dir.create(file.path(resDir, "angiosperm_subsets"))
}

### Save subset
write_csv(gbif, file.path(resDir, "angiosperm_subsets", "gbif_filtered_angiosperms_subset.csv"))

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

### Read it in again
gbif <- read_csv(file.path(resDir, "angiosperm_subsets", "gbif_filtered_angiosperms_subset.csv"))

### No. records in gbif subset Define how to split data based on number of records in each subset
records <- nrow(gbif) # 15,207,135

### Size of chunk to process each time
size <- 10000

### Define splits - 1521 chunks
splits <- seq(1, records, size)

### Number of cores to run in parallel - NOTE: on windows mcapply will only run on 1 thread
ncores <- 15

### Run parallel cleaning
mclapply(as.list(1:length(splits)),
         function(i) {
           ### Skip if already done
           outputFile <- file.path(resDir, "angiosperm_subsets", paste0("gbif_filtered_angiosperms_", i, ".csv"))

           if(!file.exists(outputFile)) {
             ### read in chunk of gbif subset
             end <- splits[i] + size - 1
             end <- ifelse((splits[i] + size - 1) > records, records, (splits[i] + size - 1))
             gbifSubset <- gbif[splits[i]:end, ]

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
subsets <- list.files(file.path(resDir, "angiosperm_subsets"), pattern = "gbif_filtered_angiosperms", full.names = TRUE)
subsets <- subsets[!grepl("gbif_filtered_angiosperms_subset.csv", subsets)]

### Make sure the column types are correct
colTypes <- read_csv(file.path(resDir, "angiosperm_subsets", "gbif_filtered_angiosperms_subset.csv"), n_max = 100000) %>%
  summarise_all(class) %>%
  as.vector() %>%
  unlist()

### Append together - 7,550,357
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
                            filterSteps      = filterStepsDup, # data frame with details of flags to apply
                            id_col           = colID,       # Column name containing unique record ID
                            lon_colname      = colLon,      # Column name containing longitude coords
                            lat_colname      = colLat,      # Column name containing latitude coords
                            pathFlaggedData  = NULL,        # path and file name of full gbif data with column indicated if flag applied. If NULL saves nothing
                            pathFilteredData = file.path(resDir, "gbif_filtered_angiosperms.csv"),  # path and file name of gbif data with flagged data filtered out.
                            plot = FALSE)

#==================================================================================================#
#--------------------------------------------- Tidy up --------------------------------------------#
#==================================================================================================#

rm(list = ls())
gc()
.rs.restartR()

