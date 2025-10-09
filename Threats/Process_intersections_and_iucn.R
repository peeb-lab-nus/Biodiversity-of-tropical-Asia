####################################################################################################
### 
### Charlie Marsh
### charliem2003@github
### 05/2025
###
### Combines the intersections and harmonises names with IUCN categories and threats
###
### To retrieve threats you will need a IUCN Red List API authentication token (https://api.iucnredlist.org).
### See the rredlist package documentation for how to add this to your r environment, which will
### depend on your OS.
###
### Taxonomic files:
### IUCN 2025:     plants, bony_fish, freshwater_fish, corals, all non-chordates
### IUCN 2024:     amphibians, freshwater_crabs, mammals, reptiles
### BirdLife 2024: birds
### POWO
###
### IUCN assessments - save in folder with structure 'taxon_name/IUCN_v2024' etc:
### IUCN 2025:     plants, bony_fish, freshwater_fish, corals, all non-chordates
### IUCN 2024:     amphibians, freshwater_crabs, mammals, reptiles
### BirdLife 2024: birds
###
### IUCN threats:
### 2025 threats data for all species
### 2024 threats data for: amphibians, birds, freshwater_crabs, mammals, reptiles
###
### You can retrieve species scores for each threat using the rredlist package which interfaces with
### the IUCN Red List API. For example, to get a list of everything, run:
###
### library(rredlist)
### codes <- rl_threats()$threats$code
### threats <- tibble()
### pb <- txtProgressBar(1, length(codes), style = 3)
### for(i in 1:length(codes)) {
###   setTxtProgressBar(pb, i)
###   print(codes[i])
###   res <- rl_threats(code = codes[i])
###   if(length(res$assessments) > 0) {
###     res <- data.frame(code = res$threat$code,
###                       name = res$threat$description$en,
###                       res$assessments)
###     threats <- bind_rows(threats, res)
###   }
### }
### write_csv(threats, file.path(threatDir, "all_threats.csv"))
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

### Load libraries
rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(rredlist)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir     <- "Threats"                                        # project dir
divDir      <- file.path("Diversity_and_endemism", "directory") # path to 'Diversity_and_endemism' folder (for non-IUCN intersections)
WCVPdir     <- file.path("WCVP", "data", "directory")           # dir with WCVP data
iucnDir     <- file.path("IUCN", "assessments", "directory")    # dir with IUCN range map intersections
birdlifeDir <- file.path("BirdLife", "range maps", "directory") # dir with BirdLife range maps (for taxonomy)

### You shouldn't need to adjust these folders
funDir    <- file.path(projDir, "Analysis_functions")  # dir with helper functions
interDir  <- file.path(projDir, "Intersections")       # dir with intersections
threatDir <- file.path(projDir, "IUCN_assessments")    # dir with threats information

#==================================================================================================#
#-------------------------------------- Read in summary data --------------------------------------#
#==================================================================================================#

### Taxon group information
taxonInfo  <- read_csv(file.path(projDir, "Global_totals.csv")) %>%
  select(Taxon, Group, Taxon_Higher, Type, Name, Colour) %>%
  mutate(Group = gsub(" NA", "", paste(Group, Taxon_Higher, sep = " "))) %>%
  mutate(GroupForPlot = gsub(" ", "\n", Group, fixed = TRUE)) %>%
  mutate(GroupForPlot = gsub("and\n", "& ", GroupForPlot, fixed = TRUE)) %>%
  mutate(TaxonForPlot = case_when(Name == "phasmids" ~ "Stick and leaf insects",
                                  Name == "sharks"   ~ "Sharks and rays",
                                  .default = Taxon))

### List of taxa
taxa <- taxonInfo$Name

### Bioregion information
regionInfo <- read_csv(file.path(projDir, "Bioregion_info.csv")) %>%
  mutate(LabelForPlot = gsub(" ", "\n", Label, fixed = TRUE))

#==================================================================================================#
#----------------------------------- Load in BTAS intersections -----------------------------------#
#==================================================================================================#

### Read in intersections and determine if native to tropical Asia
intersections <- tibble()
for(taxon in taxa) {
  group <- read_csv(file.path(interDir, paste0("Intersections_bioregions_", taxon, ".csv")),
                    show_col_types = FALSE) %>%
    mutate(Name = taxon) %>%
    relocate(Name)
  intersections <- bind_rows(intersections, group)
}

## There are some non-unique species names in different groups
intersections$UniqueName <- apply(intersections, 1, function(x) paste(x["Species"], x["Name"], sep = "_"))

### Get correct names and groups from Global_totals.csv
intersections <- intersections %>%
  left_join(select(taxonInfo, Name, Group, Name, GroupForPlot, TaxonForPlot), by = "Name") %>%
  select(-OutsideAsia, -AsiaNative) %>%
  relocate(Name, TaxonForPlot, Group, GroupForPlot, Species, UniqueName,  Genus, Family, Year, AsiaEndemic, Range_area)

#==================================================================================================#
#----------------- Load in IUCN taxonomies and just take matches with BTAS species ----------------#
#==================================================================================================#

### Read in table with POWO plant name matches for IUCN names
powo_names <- read.csv(file.path(threatDir, "plants_iucn_powo_names_table.csv"))

### Get plant IUCN taxonomy and match up with POWO look-up table
plants <- read_csv(file.path(iucnDir, "plants", "IUCN_v2025", "taxonomy.csv")) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) 

### Filter POWO-IUCN look-up table to just BTAS species - 16,055 IUCN names / 15,868 POWO names
powo_matches <- powo_names %>%
  filter(powo_name %in% intersections$Species[intersections$Group == "Vascular plants"])

### Join IUCN taxonomy and BTAS powo names
plants <- plants %>%
  left_join(powo_matches, by = "scientificName")%>%
  arrange(powo_name)

### 16,055 matches; 58,696 with no match
sum(!is.na(plants$powo_name)); sum(is.na(plants$powo_name))
plants <- plants %>%
  filter(scientificName %in% unique(scientificName[!is.na(powo_name)]))

### 167 POWO names match to multiple IUCN names (354 individual IUCN names) - remove
doubles <- unique(plants$powo_name[duplicated(plants$powo_name)])
length(doubles); length(plants$scientificName[plants$powo_name %in% doubles])
plants <- plants %>%
  filter(!powo_name %in% doubles)

### 15,701 matches between IUCN and BTAS POWO (out of 64,094 BTAS species)
nrow(plants)
plants <- plants %>%
  select(internalTaxonId, scientificName, powo_name, kingdomName, phylumName, className, orderName, familyName, genusName)

### For fish and corals we are using the 2025 range map taxonomies - remove duplicates (both marine and freshwater sp.)
fish <- bind_rows(read_csv(file.path(iucnDir, "bony_fish", "IUCN_v2025", "taxonomy.csv")),
                  read_csv(file.path(iucnDir, "freshwater_fish", "IUCN_v2025", "taxonomy.csv")),
                  read_csv(file.path(iucnDir, "corals", "IUCN_v2025", "taxonomy.csv"))) %>%
  select(internalTaxonId, scientificName, kingdomName, phylumName, className, orderName, familyName, genusName) %>%
  distinct()

### For invertebrates (excl. frehswater crabs and corals)
inverts <- read_csv(file.path(iucnDir, "non-chordates", "IUCN_v2025", "taxonomy.csv")) %>%
  filter(!familyName %in% c("GECARCINUCIDAE", "POTAMIDAE", "POTAMONAUTIDAE", "PSEUDOTHELPHUSIDAE", "TRICHODACTYLIDAE")) %>%
  filter(!phylumName == "CNIDARIA") %>%
  select(internalTaxonId, scientificName, kingdomName, phylumName, className, orderName, familyName, genusName)

### For some taxa we have 2024-1 rangemaps, in which case use those assessments
verts <- bind_rows(read_csv(file.path(iucnDir, "amphibians", "IUCN_v2024", "taxonomy.csv")),
                   read_csv(file.path(iucnDir, "freshwater_crabs", "IUCN_v2024", "taxonomy.csv")),
                   read_csv(file.path(iucnDir, "mammals", "IUCN_v2024", "taxonomy.csv")),
                   read_csv(file.path(iucnDir, "reptiles", "IUCN_v2024", "taxonomy.csv"))) %>%
  select(internalTaxonId, scientificName, kingdomName, phylumName, className, orderName, familyName, genusName)

### For birds, the taxonomy is the BirdLife International (note, they don't provide the internalTaxonId)
birds <- read_csv(file.path(birdlifeDir, "HBW_BirdLife List of Birds v.81.csv")) %>%
  dplyr::rename(scientificName = "Scientific name",
                familyName     = "Family name",
                orderName      = "Order") %>%
  mutate(kingdomName    = "ANIMALIA",
         phylumName     = "CHORDATA",
         className      = "AVES",
         genusName      = sapply(scientificName, function(x) strsplit(x, " ")[[1]][1])) %>%
  select(scientificName, kingdomName, phylumName, className, orderName, familyName, genusName)

### Get the internalTaxonId from the assessments (only for asian subset of species)
birdIds <- read_csv(file.path(iucnDir, "birds", "IUCN_v2024", "assessments.csv")) %>%
  select(internalTaxonId, scientificName)
birds <- birds %>%
  left_join(birdIds, by = "scientificName")

### Put them altogether
taxonomy <- bind_rows(plants, fish, inverts, verts, birds) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  mutate(phylumName = str_to_title(phylumName),
         className  = str_to_title(className),
         orderName  = str_to_title(orderName),
         familyName = str_to_title(familyName)) %>%
  select(internalTaxonId, scientificName, powo_name, kingdomName, phylumName, className, orderName, familyName, genusName)

any(is.na(taxonomy$scientificName))
any(duplicated(taxonomy$scientificName))

#==================================================================================================#
#------------------------------------ Load in IUCN assessments ------------------------------------#
#==================================================================================================#

### For fish and corals we used 2025 assessments
fish <- bind_rows(read_csv(file.path(iucnDir, "bony_fish", "IUCN_v2025", "assessments.csv")),
                  read_csv(file.path(iucnDir, "freshwater_fish", "IUCN_v2025", "assessments.csv")),
                  read_csv(file.path(iucnDir, "corals", "IUCN_v2025", "assessments.csv"))) %>%
  select(internalTaxonId, redlistCategory, redlistCriteria, yearPublished, assessmentDate) %>%
  left_join(taxonomy, by = "internalTaxonId")
any(is.na(fish$className))

### Also get IUCN 2025 assessments for invert groups
assess2025 <- read_csv(file.path(iucnDir, "non-chordates", "IUCN_v2025", "assessments.csv")) %>%
  select(internalTaxonId, redlistCategory, redlistCriteria, yearPublished, assessmentDate) %>%
  left_join(taxonomy, by = "internalTaxonId") %>%
  filter(!familyName %in% c("Gecarcinucidae", "Potamidae", "Potamonautidae", "Pseudothelphusidae", "Trichodactylidae")) %>%
  filter(!phylumName == "Cnidaria")
any(is.na(assess2025$className))

### For plants, swap IUCN name for POWO name match
plants <- read_csv(file.path(iucnDir, "plants", "IUCN_v2025", "assessments.csv")) %>%
  select(internalTaxonId, redlistCategory, redlistCriteria, yearPublished, assessmentDate) %>%
  left_join(taxonomy, by = "internalTaxonId") %>%
  filter(!is.na(powo_name)) %>%
  mutate(scientificName = powo_name) %>%
  select(-powo_name)
any(is.na(plants$className))

### For taxa with 2024-1 assessments, we don't have internalTaxonIds for the taxonomy - match by name
assess2024 <- bind_rows(read_csv(file.path(iucnDir, "amphibians", "IUCN_v2024", "assessments.csv")),
                        read_csv(file.path(iucnDir, "freshwater_crabs", "IUCN_v2024", "assessments.csv")),
                        read_csv(file.path(iucnDir, "mammals", "IUCN_v2024", "assessments.csv")),
                        read_csv(file.path(iucnDir, "reptiles", "IUCN_v2024", "assessments.csv"))) %>%
  select(scientificName, redlistCategory, redlistCriteria, yearPublished, assessmentDate) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  left_join(taxonomy, by = "scientificName") 
any(is.na(assess2024$className))

### For birds, we don't have global assessment data
birds <- read_csv(file.path(iucnDir, "birds", "IUCN_v2024", "assessments.csv")) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  ### Some name differences between IUCN and BirdLife taxonomy
  mutate(scientificName = case_when(scientificName == "Athene_granti"         ~ "Ninox_granti",
                                    scientificName == "Athene_jacquinoti"     ~ "Ninox_jacquinoti",
                                    scientificName == "Athene_malaitae"       ~ "Ninox_malaitae",
                                    scientificName == "Athene_roseoaxillaris" ~ "Ninox_roseoaxillaris",
                                    scientificName == "Strigops_habroptilus"  ~ "Strigops_habroptila",
                                    .default = scientificName)) %>%
  select(internalTaxonId, scientificName, redlistCategory, redlistCriteria, yearPublished, assessmentDate) %>%
  left_join(select(taxonomy, -internalTaxonId), by = "scientificName") 
any(is.na(birds$className))

### Put them together
iucn <- bind_rows(assess2024, plants, fish, assess2025, birds) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  rename(sis_taxon_id = internalTaxonId) %>%
  mutate(redlistCategory = case_match(redlistCategory,
                                      c("Critically Endangered", "CR")                  ~ "CR",
                                      c("Endangered", "EN")                             ~ "EN",
                                      c("Vulnerable", "VU")                             ~ "VU",
                                      c("Least Concern", "Near Threatened", "NT", "LC",
                                        "Lower Risk/near threatened",
                                        "Lower Risk/conservation dependent",
                                        "Lower Risk/least concern")                     ~ "NT",
                                      c("Data Deficient", "DD")                         ~ "DD",                   
                                      c("Extinct", "Extinct in the Wild", "EX", "EW")   ~ "EX")) %>%
  mutate(assessmentDate = year(assessmentDate)) %>%
  distinct()

any(is.na(iucn$scientificName))
any(is.na(iucn$sis_taxon_id))
any(is.na(iucn$yearPublished))
any(duplicated(iucn$scientificName))
any(duplicated(iucn$sis_taxon_id))

#==================================================================================================#
#-------------------------------------- Load in IUCN threats --------------------------------------#
#==================================================================================================#

### All 2025 threats data - 100,447 species
threatsAll <- read_csv(file.path(threatDir, "all_threats.csv"))

filter(iucn, (!sis_taxon_id %in% threatsAll$sis_taxon_id) & !scientificName %in% gsub(" ", "_", threatsAll$taxon_scientific_name))

### Try joining by sis_taxon_id - 57,174 species
threatsSis <- threatsAll %>%
  left_join(iucn, by = "sis_taxon_id") %>%
  filter(!is.na(scientificName))

### For those that don't match, join by name - 762 (sis_taxon_id) or 758 (internalTaxonId) species
threatsName <- threatsAll %>%
  filter(!sis_taxon_id %in% threatsSis$sis_taxon_id) %>%
  mutate(scientificName = gsub(" ", "_", taxon_scientific_name)) %>%
  left_join(iucn, by = "scientificName") %>%
  filter(!is.na(sis_taxon_id.y))

### Remove groups where we are using 2024 assessments - 39,164 species remaining
threats2025 <- threatsSis %>%
  filter(className   !=     "Amphibia") %>%
  filter(className   !=     "Aves") %>%
  filter(!familyName %in% c("Gecarcinucidae", "Potamidae", "Potamonautidae", "Pseudothelphusidae", "Trichodactylidae")) %>%
  filter(className   !=     "Mammalia") %>%
  filter(className   !=     "Reptilia") %>%
  select(-className, -familyName)

### Just for groups we are using 2025 assessments for
threats2025 <- threats2025 %>%
  filter(latest == TRUE) %>%
  pivot_wider(id_cols = c(sis_taxon_id, taxon_scientific_name),
              names_from = code,
              values_from = year_published,
              values_fn = ~ paste(unique(.x), collapse = "|"))

### Groups that we're using 2024 data for - 15,597 species
threats2024 <- bind_rows(read_csv(file.path(threatDir, "amphibians_threats.csv")),
                         read_csv(file.path(threatDir, "freshwater_crabs_threats.csv")),
                         read_csv(file.path(threatDir, "mammals_threats.csv")),
                         read_csv(file.path(threatDir, "reptiles_threats.csv"))) %>%
  mutate(code = gsub("[.]", "_", code)) %>%
  left_join(distinct(select(threatsAll, assessment_id, year_published)), # We need to get year published
            join_by(assessmentId == assessment_id)) %>%
  pivot_wider(id_cols = c(internalTaxonId, scientificName),
              names_from = code,
              values_from = year_published,
              values_fn = ~ paste(unique(.x), collapse = "|")) %>%
  rename(sis_taxon_id = internalTaxonId) %>%
  rename(taxon_scientific_name = scientificName)

### For birds, the BirdLife internalTaxonId does not seem to match up with the IUCN sis_taxon_id -3,804 species
### Use the scientificName to get correct sis_taxon_id and year_published from IUCN first
threatsBirds <- read_csv(file.path(threatDir, "birds_threats.csv")) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  left_join(select(iucn, scientificName, yearPublished), by = "scientificName") %>%
  mutate(year_published = case_when(is.na(yearPublished) ~ 1, .default = yearPublished)) %>%
  mutate(code = gsub("[.]", "_", code)) %>%
  pivot_wider(id_cols = c(internalTaxonId, scientificName),
              names_from = code,
              values_from = year_published,
              values_fn = ~ paste(unique(.x), collapse = "|")) %>%
  rename(sis_taxon_id = internalTaxonId) %>%
  rename(taxon_scientific_name = scientificName)

## Bind together
threats <- bind_rows(threats2025, threats2024, threatsBirds) %>%
  left_join(iucn, by = "sis_taxon_id") %>%
  relocate(sis_taxon_id, scientificName, powo_name, kingdomName, phylumName, className, orderName, familyName, genusName,
           redlistCategory, redlistCriteria, yearPublished, assessmentDate)

any(is.na(threats$sis_taxon_id))
any(is.na(threats$taxon_scientific_name))

### Species with no matches
filter(threats, is.na(scientificName))

### Any species with missing name matches should be birds, or a few european species from other groups
### The only exception is the Ke Go White-toothed Shrew (Crocidura kegoensis). Add manually
threats$scientificName[threats$taxon_scientific_name == "Crocidura kegoensis"] <- "Crocidura_kegoensis"
threats$kingdomName   [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "ANIMALIA"
threats$phylumName    [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "CHORDATA"
threats$className     [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "MAMMALIA"
threats$orderName     [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "EULIPOTYPHLA"
threats$familyName    [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "SORICIDAE"
threats$genusName     [threats$taxon_scientific_name == "Crocidura kegoensis"] <- "Crocidura"
threats$"5_3_5"       [threats$taxon_scientific_name == "Crocidura kegoensis"] <- 2024
threats$"11_1"        [threats$taxon_scientific_name == "Crocidura kegoensis"] <- 2024
threats$"12_1"        [threats$taxon_scientific_name == "Crocidura kegoensis"] <- 2024

### Remove others
threats <- filter(threats, !is.na(scientificName))

#==================================================================================================#
#----------------------------- Join IUCN assessments with threats data ----------------------------#
#==================================================================================================#

### For plants, replace scientificName with POWO name
iucn <- iucn %>%
  mutate(scientificName = case_when(!is.na(powo_name) ~ powo_name, .default = scientificName)) %>%
  select(-powo_name)

threats <- threats %>%
  mutate(scientificName = case_when(!is.na(powo_name) ~ powo_name, .default = scientificName)) %>%
  select(-powo_name, -sis_taxon_id, -kingdomName, -phylumName, -className, -orderName, -familyName, -genusName,
         -redlistCategory, -redlistCriteria, -yearPublished, -assessmentDate, -taxon_scientific_name)

### Join up the threats with IUCN assessments - 102,179 species
allIucn <- full_join(iucn, threats, by = "scientificName")

### Keep only those names that feature in our intersections - 33,968 species
allIucn <- allIucn %>%
  filter(scientificName %in% intersections$Species) %>%
  relocate(scientificName, sis_taxon_id) %>%
  rename(kingdomIUCN = kingdomName,
         phylumIUCN  = phylumName,
         classIUCN   = className, 
         orderIUCN   = orderName,
         familyIUCN  = familyName,
         genusIUCN   = genusName)

any(is.na(allIucn$sis_taxon_id))
any(is.na(allIucn$scientificName))
any(is.na(allIucn$redlistCategory))

#==================================================================================================#
#------------------------------------ Join everything together ------------------------------------#
#==================================================================================================#

### One case of species name in diptera and angiosperms, and a few species in both marine and freshwater fish
### Using combo of species and group to retain all cases
final <- left_join(intersections,
                   allIucn,
                   by = join_by(Species == scientificName)) 
table(final$TaxonForPlot, final$redlistCategory, useNA = "always")

### We have a few birds and a crab missing assessments. Check for ones
filter(final, is.na(redlistCategory))
missing <- filter(final, is.na(redlistCategory)) %>%
  select(Species, Name, TaxonForPlot) %>%
  mutate(redlistCategory = NA) %>%
  mutate(yearPublished   = NA) %>%
  filter(TaxonForPlot %in% c("Birds", "Freshwater crabs"))
table(missing$TaxonForPlot)

pb <- txtProgressBar(1, nrow(missing), style = 3)
for(i in 1:nrow(missing)) {
  setTxtProgressBar(pb, i)
  rl <- NULL
  sp <- strsplit(missing$Species[i], "_")[[1]]
  try(rl <- rredlist::rl_species(genus = sp[1], species = sp[2]), silent = TRUE)
  if(!is.null(rl)) {
    rl <- as.data.frame(rl$assessments)
    if(nrow(rl) > 0) {
      if(missing$TaxonForPlot[i] %in% c("Amphibians", "Birds", "Freshwater crabs", "Mammals", "Reptiles")) {
        rl <- rl[rl$year_published <= 2024 & grepl("Global", rl$scopes), ]
      }
      if(!missing$TaxonForPlot[i] %in% c("Amphibians", "Birds", "Freshwater crabs", "Mammals", "Reptiles")) {
        rl <- rl[rl$year_published <= 2025 & grepl("Global", rl$scopes), ]
      }
      if(nrow(rl) > 0) {
        final$redlistCategory[final$Species == missing$Species[i]] <- rl$red_list_category_code[which.max(rl$year_published)]
        final$yearPublished[final$Species == missing$Species[i]]   <- rl$year_published[which.max(rl$year_published)]
      }
    }
  }
}
table(final$TaxonForPlot, final$redlistCategory, useNA = "always")

### Also deal with 2 cases where plant species and diptera species have same name - the dipteras should be NA
final[final$Species %in% c("Limnophila_glabra", "Ormosia_formosana") & final$Name == "diptera",
      which(names(final) == "sis_taxon_id"):ncol(final)] <- NA

### Final tidying up
allThreats <- final %>%
  mutate(redlistCategory = case_when(is.na(redlistCategory)  ~ "Not evaluated",
                                     redlistCategory == "LC" ~ "Non-threatened",
                                     redlistCategory == "NT" ~ "Non-threatened",
                                     .default = redlistCategory)) %>%
  mutate(yearPublished = as.numeric(yearPublished)) %>%
  mutate(TaxonForPlot = factor(TaxonForPlot, levels = unique(taxonInfo$TaxonForPlot))) %>%
  mutate(GroupForPlot = case_when(
    Group == "Vascular plants" ~ "\nVascular\nplants",
    Group == "Terrestrial and freshwater vertebrates" ~ "Terrestrial\n& freshwater\nvertebrates",
    Group == "Terrestrial and freshwater invertebrates" &
      Name %in% c("ants", "bees", "butterflies", "trichoptera", "chilopods",
                  "diptera") ~ "Terrestrial\n& freshwater\ninvertebrates",
    Group == "Terrestrial and freshwater invertebrates" &
      Name %in% c("freshwater_crabs", "diplopods", "miridae", "spiders",
                  "phasmids") ~ "Terrestrial\n& freshwater\ninvertebrates",
    Group == "Marine" ~ "\n\nMarine")) %>%
  mutate(GroupForPlot = factor(GroupForPlot,
                               levels = c("\nVascular\nplants",
                                          "Terrestrial\n& freshwater\nvertebrates",
                                          "Terrestrial\n& freshwater\ninvertebrates",
                                          "\n\nMarine")
  ))

table(allThreats$Name, allThreats$redlistCategory, useNA = "always")
apply(allThreats[, 1:45], 2, function(x) sum(is.na(x)))

### No. sp with RL categories but no threats listed
threatData <- allThreats %>%
  mutate(across(`1_1`:`12_1`, ~ !is.na(.))) %>%
  mutate(threatData = rowSums(.[,  which(names(allThreats) == "1_1"):which(names(allThreats) == "12_1")])) %>%
  mutate(threatData = case_when(threatData > 0 ~ 1, .default = 0)) %>%
  relocate(threatData)
table(threatData$redlistCategory, threatData$threatData)

### Save final threats info
write_csv(allThreats, file.path("All_iucn_and_threats.csv"))
