####################################################################################################
### Intersections between Miridae and bioregions
### Kendrick Chan Siu Chun and Charlie Marsh
### charliem2003@github
### 07/2024
###
### Taxonomy used:    Planetary Biodiversity Inventory (PBI) for Plant Bugs
### Checklist used:   Planetary Biodiversity Inventory (PBI) for Plant Bugs
###
### CITATION:
### Schuh, R. T. (2016). On-line Systematic Catalog of Plant Bugs (Insecta - Heteroptera - Miridae)
###  (Mar 2013). http://research.amnh.org/pbi/index.html
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

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity"                                  # project dir
dataDir   <- file.path("PBI", "data", "directory")        # dir that contains the checklist data

### You shouldn't need to adjust these folders
resDir    <- file.path(projDir, "Intersections")          # dir to save results to
funDir    <- file.path(projDir, "Analysis_functions")     # dir that contains the function scripts
lookupDir <- file.path(resDir, "Look-up_tables")          # dir that contains look-up table for assigning bioregions

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for cleaning dates and extracting region info from inside brackets
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### CoL Hemipterans - Miridae Data
col_miridae <- read.csv(file.path(dataDir, "CoL_Hemiptera.csv"))%>%
  as_tibble()

### Spreadsheet with codes to assign to countries and provinces to bioregions
bioregions <- read.csv(file.path(lookupDir, "country_codes_BTAS_miridae.csv")) %>%
  select(Regions, Bioregion, AsiaNative)

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### Make column for species name and region info - change for relevant column names
col_region <- col_miridae %>%
  mutate(Species = gsub(" ", "_", col.scientificName)) %>%
  mutate(Regions = tolower(col.area)) %>%
  rename(Genus = col.genericName) %>%
  mutate(Family = "Miridae") %>%
  mutate(Year = clean_publication_dates(col.authorship)) %>%
  select(Species, Genus, Family, Regions, Year)

### Separate out individual entries to unique rows
col_region <- col_region %>%
  separate_longer_delim(Regions, delim = "; ")

### Regions within counties are in the form Country: province1, province2
### This is a  version of the extractParentheses function to separate them out to individual entries
extractColon <- function(distribution) {
  splits <- strsplit(distribution, ": |, ")[[1]]
  splits <- as.vector(sapply(splits, function(x) trimws(x, "both")))
  if(length(splits) > 1) {
    splits <- paste(splits[1], splits[2:length(splits)], sep = ": ")
  }
  final <- paste(splits, collapse = ", ")
  return(final)
}

### Apply function and put on individual rows
col_region <- col_region %>%
  mutate(Regions = sapply(Regions, function(x) extractColon(x))) %>%
  separate_longer_delim(Regions, delim = ", ")

### Out of 836 species with records in China and/or Japan there are 68 species in raw data with
### only 'China' or 'Japan' and no other province info -> retain China and Japan info
CN_JP_provinces <- c("amami", "china:", "xinjiang", "inner mongolia", "fujian", "gansu", "guangdong",
                     "guangxi", "hainan", "hebei", "heilongjiang", "henan", "hokkaido", "hong kong",
                     "hongkong", "honshu", "hubei", "jiangxi", "kunashir", "kyushu", "manchuria",
                     "ningxia", "north china", "ogasawara", "okinawa", "rishiri", "ryuku", "ryukyu",
                     "s china", "saipan", "shaanxi", "shandong", "shanxi", "shikoku", "sichuan",
                     "tibet", "xinjiang", "xizang", "yunnan", "zhejiang")
length(unique(col_region$Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), col_region$Regions)]))
col_region %>%
  filter(Species %in% col_region$Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), col_region$Regions)]) %>%
  filter(!Species %in% col_region$Species[grepl(paste(CN_JP_provinces, collapse = "|"), col_region$Regions)]) %>%
  summarise(length(unique(Species)))

### If country == Indonesia or Malaysia, but more detailed province info also given, remove those rows
spMalaysia <- unique(col_region$Species[col_region$Species %in% col_region$Species[col_region$Regions %in% c("malaysia", "malaya")]])
for(i in 1:length(spMalaysia)) {
  sp <- filter(col_region, Species == spMalaysia[i])
  if(sum(grepl(paste(c("malaysia:", "malaya:", "malacca", "penang", "west malaysia", "boreno",
                       "borneo", "bruei", "brunei", "east malaysia", "sabah", "pulo laut", "sarawak"),
                     collapse = "|"), sp$Regions)) >= 1) {
    col_region <- filter(col_region, Species != spMalaysia[i] | !Regions %in% c("malaysia", "malaya"))
  }
}

spIndonesia <- unique(col_region$Species[col_region$Species %in% col_region$Species[col_region$Regions == "indonesia"]])
for(i in 1:length(spIndonesia)) {
  sp <- filter(col_region, Species == spIndonesia[i])
  if(sum(grepl(paste(c("indonesia:", "boreno", "borneo", "bruei", "brunei", "pulo laut", "sarawak",
                       "java", "lombok", "timor", "amboina", "ambon", "moluccas", "ternate", "bismark",
                       "madang", "new britain", "new guinea", "new ireland", "new guinea", "sulawesi",
                       "batu", "engano", "mantawei", "si-oban", "sumatra"),
                     collapse = "|"), sp$Regions)) >= 1) {
    col_region <- filter(col_region, Species != spIndonesia[i] | !Regions == "indonesia")
  }
}

### assign to BTAS bioregions
col_region <- col_region %>%
  left_join(bioregions, by = "Regions", relationship = "many-to-many") %>%
  mutate(AsiaNative = case_when(AsiaNative == 0   ~ 0,
                                AsiaNative == 1   ~ 1,
                                is.na(AsiaNative) ~ 0))

### which species are native to tropical asia
native <- col_region %>%
  group_by(Species) %>%
  reframe(native = max(AsiaNative, na.rm = TRUE)) %>%
  mutate(native = case_when(native >= 1 ~ 1,
                            native == 0 ~ 0))
col_region <- col_region %>%
  left_join(native, by = "Species") %>%
  filter(native == 1) %>%
  distinct() %>%
  mutate(Bioregion = case_when(is.na(Bioregion) & AsiaNative == 1 ~ "Unknown",
                               .default = Bioregion))

### which species are endemic to tropical asia
endemic <- col_region %>%
  group_by(Species) %>%
  reframe(AsiaEndemic = any(is.na(Bioregion))) %>%
  mutate(AsiaEndemic = case_match(AsiaEndemic,
                                  FALSE ~ 1,
                                  TRUE  ~ 0))
col_region <- left_join(col_region, endemic, by = "Species")

### Filter only native intersections
col_region <- col_region %>%
  filter(AsiaNative == 1)

### convert to species x bioregion data frame
col_region <- col_region %>%
  select(Species, Genus, Family, Year, Bioregion, AsiaNative, AsiaEndemic) %>%
  filter(!duplicated(.))  %>%
  pivot_wider(id_cols = c(Species, Genus, Family, Year, AsiaEndemic),    
              names_from  = Bioregion,
              values_from = AsiaNative, 
              values_fill = 0) %>%
  mutate(AsiaNative = 1) %>%
  relocate(Species, Genus, Family, Year, AsiaNative, AsiaEndemic)

### We have 12 species where tropical asia distribution is unknown
sum(col_region$Unknown)

### But there are 11 species where bioregion is Unknown but the species also has other bioregion info
### In these cases we will remove the Unknown definition
bioregions <- unique(bioregions$Bioregion)
bioregions <- bioregions[!is.na(bioregions)]
col_region$Unknown[col_region$Unknown == 1 & rowSums(col_region[, bioregions]) > 0] <- 0

### total number of species native to region and number of species endemic to tropical asia
c(native = nrow(col_region),
  endemic = sum(col_region$AsiaEndemic))
colSums(col_region[, c(bioregions, "Unknown")])

### save results
write.csv(col_region, file.path(resDir, "Intersections", "Intersections_bioregions_miridae.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
