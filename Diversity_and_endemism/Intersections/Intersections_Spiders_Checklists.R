####################################################################################################
### Run intersections between World Spider Catalog (WSC) checklists and bioregions
### Charlie Marsh and original code from Anna Holmquist
### charliem2003@github
### 05/2024
###
### Taxonomy used:    World Spider Catalog
### Checklist used:   World Spider Catalog
###
### CITATION: World Spider Catalog (2024). World Spider Catalog. Version 25.0. Natural History Museum Bern,
###  online at http://wsc.nmbe.ch, accessed on 19/02/2024. doi: 10.24436/2
###
### saves csv with the presence-absence of each species within each bioregion
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(dplyr)
library(tidyr)

### Locations of data, scripts and results
baseDir   <- "/mnt/Work"
regionDir <- file.path(baseDir, "NUS", "BTAS_data", "Bioregions") # dir that contains the regions data
dataDir   <- file.path(baseDir, "spatial_data", "biodiversity", "checklists", "spiders") # dir that contains the checklist data
resDir    <- file.path(baseDir, "NUS", "BTAS", "Intersections") # dir to save results to
funDir    <- file.path(baseDir, "NUS", "BTAS", "Analysis_functions")   # dir that contains the function scripts


### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity"                                  # project dir
dataDir   <- file.path("MilliBase", "data", "directory")  # dir that contains the checklist data

### You shouldn't need to adjust these folders
resDir    <- file.path(projDir, "Intersections")          # dir to save results to
funDir    <- file.path(projDir, "Analysis_functions")     # dir that contains the function scripts
lookupDir <- file.path(resDir, "Look-up_tables")          # dir that contains look-up table for assigning bioregions

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for extracting region info from inside brackets
source(file.path(funDir, "Intersections", "extractParantheses.R"))

### World Spider Catalog Data
wsc <- read.csv(file.path(dataDir, "World_Spider_Catalog", "species_export_20240313.csv")) %>%# "WSC_20240219.csv")) %>%
  as_tibble() %>% 
  mutate(Species = paste(genus, species, sep = "_"))

### There is one subspecies (not in region) and a species can have multiple rows, so n species = 51528
wsc$Species[grepl(" ", wsc$Species)]

### Spreadsheet with codes to assign to countries and provinces to bioregions
bioregions <- read.csv(file.path(lookupDir, "country_codes_BTAS_spiders.csv")) %>%
  select(Region, Bioregion, AsiaNative)

### Spreadsheet outlining bioregions when distribution information is 'X to Y'
paths <- read.csv(file.path(lookupDir, "country_paths_BTAS_spiders.csv")) %>%
  select(Raw, Region, AsiaNative) %>%
  filter(!is.na(Region)) %>%
  separate_longer_delim(Region, delim = ", ") %>%
  rename(Bioregion = Region) %>%
  rename(Region = Raw)

### Join together to get full list of names
bioregions <- bind_rows(bioregions, paths)

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### Remove Introduced countries (assuming always comes at list at end) and questions
# wsc$distribution[grepl("Introduced", wsc$distribution)]
wsc_region <- wsc %>% 
  mutate(Region = sapply(distribution, function(x) strsplit(x, ". Introduced")[[1]][1])) %>%
  filter(grepl("[?|possibly]", Region))

### Extract province/island info from brackets and append to country name
wsc_region <- wsc_region %>%
  mutate(Region = sapply(Region, extractParantheses))

### Separate out individual entries to unique rows
wsc_region <- wsc_region %>%
  separate_longer_delim(Region, delim = ", ")

### Out of 7115 species with records in China and/or Japan there are 6639 species in raw data with
### only 'China' or 'Japan' and no other province info -> discard all China and Japan info
CN_JP <- c("China", "to China", "China to", "and China", "China and", "Cina",
           "Japan", "to Japan", "Japan to", "and Japan", "Japan and")
CN_JP_provinces <- c("Hainan", "China Himalaya", "Hong Kong", "Inner Mongolia", "Tibet", "Yarkand",
                     "Yongxing", "Himalaya China", "Iriomote", "Japan mainland", "Ogasawara", "Okinawa",
                     "Ryukyu", "Tsushima")
length(unique(wsc_region$Species[grepl(paste(CN_JP, collapse = "|"), wsc_region$distribution)]))
wsc_region %>%
  filter(Species %in% wsc_region$Species[grepl(paste(c(CN_JP, CN_JP_provinces), collapse = "|"), wsc_region$distribution)]) %>%
  filter(!Species %in% wsc_region$Species[grepl(paste(CN_JP_provinces, collapse = "|"), wsc_region$distribution)]) %>%
  summarise(length(unique(Species)))

### If country == Indonesia or Malaysia, but more detailed province info also given, remove those rows
spMalaysia <- unique(wsc_region$Species[wsc_region$Species %in% wsc_region$Species[wsc_region$Region == "Malaysia"]])
for(i in 1:length(spMalaysia)) {
  sp <- filter(wsc_region, Species == spMalaysia[i])
  if(sum(grepl(paste(c("Borneo", "Brunei", "Gaya", "Malaysia mainland", "peninsula", "Peninsula",
                       "peninsular", "penisula", "Tioman"),
                     collapse = "|"), sp$Region)) >= 1) {
    wsc_region <- filter(wsc_region, Species != spMalaysia[i] | !Region == "Malaysia")
  }
}

spIndonesia <- unique(wsc_region$Species[wsc_region$Species %in% wsc_region$Species[wsc_region$Region == "Indonesia"]])
for(i in 1:length(spIndonesia)) {
  sp <- filter(wsc_region, Species == spIndonesia[i])
  if(sum(grepl(paste(c("Borneo", "Gaya", "Kalimantan", "Sunda", "Bali", "Java", "Flores", "Lombok",
                       "Sumba", "Sumbawa", "Timor", "Amboina", "Ambon", "Aru", "Banda", "Buru", "Halmahera",
                       "Kei", "Moluccas", "Maluku", "Seram", "Ternate", "Bismarck", "Biak", "New Guinea",
                       "West Papua", "New Ireland", "New Britain", "Sulawesi", "Togian", "Sumatra",
                       "Belitung", "Bodjo", "Krakatau", "Krakatoa", "Mentawai", "Nias", "Simeulue"),
                     collapse = "|"), sp$Region)) >= 1) {
    wsc_region <- filter(wsc_region, Species != spIndonesia[i] | !Region == "Indonesia")
  }
}

### assign to BTAS bioregions
wsc_region <- wsc_region %>%
  left_join(bioregions, by = "Region", relationship = "many-to-many") %>%
  select(family, genus, year, distribution, Region, Species, Bioregion, AsiaNative)

### to check which countries aren't matched to Bioregions
# sort(unique(select(filter(wsc_region, is.na(Bioregion)), Region)$Region))

## which species are native to tropical asia
native <- wsc_region %>%
  group_by(Species) %>%
  reframe(native = max(AsiaNative, na.rm = TRUE)) %>%
  mutate(native = case_when(native >= 1 ~ 1,
                            native == 0 ~ 0))
wsc_region <- wsc_region %>%
  left_join(native, by = "Species") %>%
  filter(native == 1) %>%
  distinct() %>%
  mutate(Bioregion = case_when(is.na(Bioregion) & AsiaNative == 1 ~ "Unknown",
                               .default = Bioregion))

### which species are endemic to tropical asia
endemic <- wsc_region %>%
  group_by(Species) %>%
  reframe(AsiaEndemic = any(is.na(Bioregion))) %>%
  mutate(AsiaEndemic = case_match(AsiaEndemic,
                                  FALSE ~ 1,
                                  TRUE  ~ 0))
wsc_region <- left_join(wsc_region, endemic, by = "Species")

### Filter only native intersections
wsc_region <- wsc_region %>%
  filter(AsiaNative == 1)

### convert to species x bioregion data frame
wsc_region <- wsc_region %>%
  select(Species, genus, family, year, Bioregion, AsiaNative, AsiaEndemic) %>%
  filter(!duplicated(.))  %>%
  pivot_wider(id_cols = c(Species, genus, family, year, AsiaEndemic),    
              names_from  = Bioregion,
              values_from = AsiaNative, 
              values_fill = 0) %>%
  mutate(AsiaNative = 1) %>%
  dplyr::rename(Genus = genus, Family = family, Year = year) %>%
  relocate(Species, Genus, Family, Year, AsiaNative, AsiaEndemic)

### Sometimes subspecies are given a unique speciesId so a species may have multiple rows. Merge them
bioregions <- unique(bioregions$Bioregion)[!is.na(unique(bioregions$Bioregion))]
wsc_region <- wsc_region %>%
  group_by(Species) %>%
  reframe(Species     = unique(Species),
          Genus       = unique(Genus),
          Family      = unique(Family),
          Year        = min(Year),
          AsiaEndemic = min(AsiaEndemic),
          across(any_of(bioregions), ~ max(.)),
          Unknown     = max(Unknown))

### We have 17 species where tropical asia distribution is unknown
sum(wsc_region$Unknown)

### But there are 5 cases where bioregion is Unknown but the species also has other bioregion info
### In these cases we will remove the Unknown definition
wsc_region$Unknown[wsc_region$Unknown == 1 & rowSums(wsc_region[, bioregions])] <- 0

### total number of species native to region and number of species endemic to tropical asia
c(native = nrow(wsc_region),
  endemic = sum(wsc_region$AsiaEndemic))
colSums(wsc_region[, c(bioregions, "Unknown")])

### save results
write.csv(wsc_region, file.path(resDir, "Intersections", "Intersections_bioregions_spiders.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
