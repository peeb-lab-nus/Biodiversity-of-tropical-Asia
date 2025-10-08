####################################################################################################
### Run intersections between John Ascher's bees checklists and bioregions
### Charlie Marsh
### charliem2003@github
### 05/2024
###
### Taxonomy used:  Ascher & Pickering 2025. Discover Life bee species guide and world checklist. Draft 58
### Checklist used: Ascher & Pickering 2025. Discover Life bee species guide and world checklist. Draft 58
###
### CITATION:
### Ascher, J. S. & Pickering, J. 2025. Discover Life bee species guide and world checklist,
###   (Hymenoptera: Apoidea: Anthophila). Draft 58. https://www.discoverlife.org/mp/20q?guide=Apoidea_species
### Orr, M. C., Hughes, A. C., Chesters, D., Pickering, J., Zhu, C.-D., & Ascher, J. S. (2021).
###   Global patterns and drivers of bee distribution. Current Biology, 31(3), 451-458.e4.
###   https://doi.org/10.1016/j.cub.2020.10.053
###
### saves csv with the presence-absence of each species within each bioregion
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
projDir   <- "Diversity_and_endemism"                     # project dir
dataDir   <- file.path("Bee", "data", "directory")        # dir that contains the checklist data

### You shouldn't need to adjust these folders
resDir    <- file.path(projDir, "Intersections")          # dir to save results to
funDir    <- file.path(projDir, "Analysis_functions")     # dir that contains the function scripts
lookupDir <- file.path(resDir, "Look-up_tables")          # dir that contains look-up table for assigning bioregions

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Function for cleaning dates and extracting region info from inside brackets
source(file.path(funDir, "Discovery_rates", "clean_publication_dates.R"))

### John's checklist data
bees <- read_csv(file.path(dataDir, "Ascher bee data tropical Asia.csv"))
colnames(bees) <- gsub(" ", "_", colnames(bees))

### Provinces code
codes <- read.csv(file.path(dataDir, "Ascher primary divisions states provinces islands codes_bioregions.csv"))

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### data frame that contains bioregion info for the bee distribution codes
codes <- codes %>%
  mutate(Code = paste(country, abv, sep = ":")) %>%
  select(Country, Bioregion, Code)

### full list of distribution codes for tropical asia
countries <- unique(codes$Code)
bioregions <- unique(codes$Bioregion)
bioregions <- bioregions[!is.na(bioregions)]

### Distribution codes in the checklist are in the following format
### country1:province1;province2 country2: country3:province
### first filter out synonyms which only have duplicated information
bees <- bees %>%
  filter(Status == "Valid species") %>%
  select(Family, Genus, Date, valid_genus_species_subspecies, Global_distribution_in_Code, SE_Asia_paper)

### loop through and convert species with multiple distribution codes to long format
spList <- unique(bees$valid_genus_species_subspecies)
intersections <- tibble()
for(i in spList) {
  sp <- filter(bees, valid_genus_species_subspecies == i)
  
  ### There are a couple of species with multiple rows
  if(nrow(sp) > 1) {
    sp <- sp %>%
      reframe(Family = unique(Family),
                  Genus  = unique(Family),
                  Date   = min(Date),
                  valid_genus_species_subspecies = unique(valid_genus_species_subspecies),
                  Global_distribution_in_Code = paste0(Global_distribution_in_Code, collapse = ";"),
                  SE_Asia_paper = case_when(any(SE_Asia_paper %in% c("Present core area",
                                                                     "Present in Pakistan",
                                                                     "Endemic to area if Pakistan included")) ~ "Present core area",
                                            all(SE_Asia_paper == "Endemic species core area") ~ "Endemic species core area"),
              .by = "valid_genus_species_subspecies")
  }
  
  ### split distribution code into countries
  spCode <- unlist(strsplit(sp$Global_distribution_in_Code, split = " ")) %>%
    strsplit(";")
  
  ### clean up some formatting/spelling errors
  spCode <- lapply(spCode, function(x) gsub("[?]", "?", x)) # only 2 cases
  spCode <- lapply(spCode, function(x) gsub("::", ":", x))
  spCode <- lapply(spCode, function(x) gsub("]", "", x))
  
  ### where multiple provinces per country append country name to make it country:province
  fullList <- vector()
  fullList <- lapply(spCode, function(x) {
    if(length(x) > 1) {
      ctry <- paste0(strsplit(x, ":")[[1]][1], ":")
      provinces <- c(x[1], paste0(ctry, x[2:length(x)]))
    } else {
      provinces <- x
    }
    fullList <- c(fullList, provinces)
    return(fullList)
  })
  fullList <- unlist(fullList)
  fullList <- gsub(" ", "", fullList)
  
  ### append to make long format data frame
  intersections <- bind_rows(intersections,
                   data.frame(Species      = sp$valid_genus_species_subspecies,
                              Genus        = sp$Genus,
                              Family       = sp$Family,
                              Date         = sp$Date,
                              AsiaEndemic  = sp$SE_Asia_paper,
                              Code         = fullList))
}

### join with code spreadsheet to assign to bioregions
intersections <- left_join(intersections, codes, by = "Code", relationship = "many-to-many")

### Out of 485 species with records in China there are only 3 species in raw data with only 'China'
### and no other china province info 
length(unique(intersections$Species[grepl("CH:", intersections$Code, fixed = FALSE)]))
intersections %>%
  filter(Species  %in% intersections$Species[intersections$Code == "CH:"]) %>%
  filter(!Species %in% intersections$Species[grepl(
    paste0(paste0("CH:", paste0(expand.grid(LETTERS, LETTERS)[, "Var1"],
                                expand.grid(LETTERS, LETTERS)[, "Var2"]), collapse = "|"),
           "|", "HK:", "|",
           paste0("HK:", paste0(expand.grid(LETTERS, LETTERS)[, "Var1"],
                                expand.grid(LETTERS, LETTERS)[, "Var2"]), collapse = "|")),
  intersections$Code)]) %>%
  summarise(length(unique(Species)))

### There are 0 species in raw data with only 'Japan'
intersections %>%
  filter(Species  %in% intersections$Species[intersections$Code == "JP:"]) %>%
  nrow()
#***# -> We're happy to keep yunnan, guangxi, okinawa etc in data

### where there is no match then it is outside our tropical asia definition
intersections <- intersections %>%
  filter(!is.na(Bioregion))

### convert to species (rows) x bioregions (cols) presence-absence table
intersections <- intersections %>%
  mutate(presence = 1) %>%
  pivot_wider(id_cols = c(Species, Genus, Family, Date, AsiaEndemic),    
              names_from = Bioregion,
              values_from = presence, 
              values_fill = 0,
              values_fn   = max) %>%
  mutate(Date = clean_publication_dates(Date, max = TRUE, first = FALSE)) %>%
  dplyr::rename(Year = Date) %>%
  mutate(AsiaEndemic = case_match(AsiaEndemic,
                                  "Endemic species core area"            ~ 1,
                                  "Present core area"                    ~ 0,
                                  "Endemic to area if Pakistan included" ~ 0))

### total number of species native to region and number of species endemic to tropical asia
c(native = nrow(intersections),
  endemic = sum(intersections$AsiaEndemic == 1))

### Save results  
write.csv(intersections, file.path(resDir, "Intersections", "Intersections_bioregions_bees.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
