####################################################################################################
### Extracts points for Porifera from OBIS species and matches to WoRMS names using 
### 
### Charlie Marsh
### charliem2003@github
### 12/2024
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(readr)
require(dplyr)
library(duckdb)
library(arrow)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir  <- "Diversity"                                                    # project dir
obisDir  <- file.path("OBIS", "data dump", "directory")                    # dir that contains OBIS data dump
wormsDir <- file.path("WoRMS", "data dump", "directory")                   # dir that contains OBIS data dump

### You shouldn't need to adjust these folders
resDir  <- file.path(projDir, "Intersections", "Intersections", "Sponges") # dir to save results to
funDir  <- file.path(projDir, "Analysis_functions")                        # dir that contains the function scripts

### Create folders for storing species intersections
if(!dir.exists(resDir)) { dir.create(resDir, recursive = TRUE) }

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for extracting description year from authorships
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### Full OBIS data for getting coordinates etc of final choice of records
con <- DBI::dbConnect(duckdb(), config = list("memory_limit"="156GB"))
obis_par <- open_dataset(file.path(obisDir, "obis_20240723.parquet"))

### Subset Porifera from matched OBIS data
obis <- to_duckdb(obis_par, con) %>%
  filter(phylum == "Porifera") %>%
  select(id, species, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, coordinatePrecision,
         depth, marine, brackish, freshwater, terrestrial, basisOfRecord, occurrenceStatus, typeStatus,
         continent, country, island, stateProvince, locality, absence,
         eventDate, year, month, day, redlist_category, dropped, absence, flags) %>%
  as_tibble()

### Read in WoRMS data for matching names and getting global richness
worms <- read_delim(file.path(wormsDir, "taxon.txt"), delim = "\t") %>%
  filter(taxonRank %in% c("Species")) %>%
  filter(taxonomicStatus == "accepted") %>%
  mutate(descriptionYear = clean_publication_dates(scientificNameAuthorship, max = TRUE)) %>%
  mutate(descriptionYear = as.numeric(trimws(descriptionYear, "both"))) %>%
  mutate(descriptionYear = case_when(descriptionYear < 1753       ~ NA,
                                     descriptionYear > 2024       ~ NA,
                                     is.infinite(descriptionYear) ~ NA,
                                     .default = descriptionYear)) %>%
  mutate(descriptionYear = case_when(is.na(scientificNameAuthorship) ~ NA,
                                     .default = descriptionYear)) %>%
  group_by(acceptedNameUsage) %>%
  filter(phylum == "Porifera") %>%
  filter(descriptionYear     == min(descriptionYear, na.rm = TRUE)) %>%
  filter(acceptedNameUsageID == min(acceptedNameUsageID, na.rm = TRUE)) %>%
  select(acceptedNameUsage, acceptedNameUsageID, order, family, genus, specificEpithet, descriptionYear)

any(duplicated(worms$acceptedNameUsage)); any(duplicated(worms$acceptedNameUsageID))

#==================================================================================================#
#------------------------------------ Match OBIS data to WoRMS ------------------------------------#
#==================================================================================================#

### Join with worms to get description year
sponges <- left_join(worms, obis,
            join_by("acceptedNameUsage" == "species"),
            copy = TRUE) %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  rename(species = acceptedNameUsage) %>%
  select(id, species, acceptedNameUsageID, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters,
         marine, brackish, freshwater, terrestrial, occurrenceStatus, absence,
         coordinatePrecision, specificEpithet, genus, family, order, descriptionYear)

### Save
write_csv(sponges, file.path(resDir, "obis_worms_sponges.csv"))

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
