####################################################################################################
### Intersections between Sponges and marine ecooregions
### Checklist generated from cleaned OBIS data by Buntaro Kusumoto
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### Taxonomy used:    WoRMS (downloaded 2024-11-11)
### Checklist used:   OBIS (occurrence.parquet - official data dump 23/07/2024)
###
### NOTE: First you must run the scripts (in 'Analysis_functions/Taxon_specific/Sponges'):
###       1) 'match_obis_to_worms_sponges.R'
###          - filters sponges from OBIS data dump, matches names to WoRMS, generates obis_worms_sponges.csv
###       2) '009_create_PresenceAbsenceTable_sponges.R'
###          - generates obis_sp_province_matrix_sponge.csv
###
### CITATION:
###   OBIS (2024) Ocean Biodiversity Information System. Intergovernmental Oceanographic Commission
###     of UNESCO. https://obis.org.
###   WoRMS Editorial Board (2024). World Register of Marine Species. Available from
###     https://www.marinespecies.org at VLIZ. Accessed 2024-11-11. doi:10.14284/170 
### 
### saves csv with the presence-absence of each species within each ecooregion
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(readr)
library(dplyr)
library(tidyr)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity_and_endemism"                       # project dir

### You shouldn't need to adjust these folders
resDir  <- file.path(projDir, "Intersections")              # dir to save results to
dataDir <- file.path(resDir,  "Intersections", "Sponges")   # dir that contains the checklist data

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Read in obis intersections
all <- read.csv(file.path(dataDir, "obis_sp_province_matrix_sponge.csv")) %>%
  as_tibble()

ecoregions <- read_csv(file.path(projDir, "Bioregion_info.csv")) %>%
  filter(Realm == "Marine") %>%
  select(Bioregion) %>%
  unlist() %>%
  as.vector()

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### The intersections have been calculated, we just need to do some tidying up
all <- all %>%
  mutate(Species = gsub(" ", "_", species)) %>%
  rename(Genus   = genus) %>%
  rename(Family  = family) %>%
  rename(Year    = descriptionYear) %>%
  rename_with(., ~ gsub(".", "_", .x, fixed = TRUE)) %>%
  select(-species, -specificEpithet, -order) %>%
  relocate(Species, Genus, Family, Year, AsiaEndemic, OutsideAsia) %>%
  distinct()

### 8,731 species globally
length(unique(all$Species))

### Filter out species not native to tropical Asia - 1,442 species
all <- all %>%
  filter(rowSums(.[, ecoregions]) > 0)

### save results
write.csv(all, file.path(resDir, "Intersections", "Intersections_bioregions_sponges.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
