####################################################################################################
### GET PLANT THREATS
### Harmonises POWO names to IUCN names
### 05/2025
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

# PACKAGES ============
rm(list = ls())
library(dplyr)
library(tidyr)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir     <- "Threats"                                        # project dir
WCVPdir     <- file.path("WCVP", "data", "directory")           # dir with WCVP data

### You shouldn't need to adjust these folders
interDir  <- file.path(projDir, "Intersections")       # dir with intersections
threatDir <- file.path(projDir, "IUCN_assessments")    # dir with threats information

#==================================================================================================#
#----------------------------------------- Read in data -------------------------------------------#
#==================================================================================================#

### Plant intersections with tropical Asia
plant_distr <- bind_rows(data.frame(Taxon = "Flowering plants",
                                    read.csv(file.path(interDir, "Intersections_bioregions_angiosperms.csv"))),
                         data.frame(Taxon = "Ferns",
                                    read.csv(file.path(interDir, "Intersections_bioregions_ferns.csv"))),
                         data.frame(Taxon = "Gymnosperms",
                                    read.csv(file.path(interDir, "Intersections_bioregions_gymnosperms.csv"))),
                         data.frame(Taxon = "Lycophytes",
                                    read.csv(file.path(interDir, "Intersections_bioregions_lycophytes.csv"))))

### IUCN threats
plant_threats <- read.csv(file.path(threatDir, "all_threats.csv")) %>%
  select(code, name, year_published, latest, sis_taxon_id) %>%
  mutate(presence = 1) %>%
  distinct() %>%
  pivot_wider(id_cols = c(year_published, latest, sis_taxon_id),
              names_from = code,
              values_from = presence)

### IUCN assessments
plant_assessments <- read.csv(file.path(iucnDir, "plants", "IUCN_v2025", "assessments.csv")) %>%
  select(internalTaxonId, scientificName, redlistCategory, redlistCriteria, yearPublished,
         assessmentDate, systems, realm, yearLastSeen, possiblyExtinct, possiblyExtinctInTheWild, scopes)
sum(duplicated(plant_assessments$scientificName)) # no duplicates
length(unique(plant_assessments$scientificName))  # 74751 species globally

### Add threats to assessments
plant_assessments <- left_join(plant_assessments,
                               plant_threats,
                               join_by("internalTaxonId" == "sis_taxon_id"))
filter(plant_assessments, is.na(year_published)) # 7 species have no threats info

### POWO names
plant_names <- read.table(file.path(WCVPdir, "wcvp_names.csv"),
                          header = TRUE, sep = "|", quote = "", fill = TRUE, encoding = "UTF-8")
sum(plant_assessments$scientificName %in% plant_names$taxon_name) # 73789 species
table(subset(plant_names, taxon_name %in% plant_assessments$scientificName)$taxon_status) # Only 68777 accepted

#==================================================================================================#
#------------------------------- Reconcile IUCN names against POWO --------------------------------#
#==================================================================================================#

plant_assessment_names <- tibble(scientificName = plant_assessments$scientificName) %>%
  left_join(select(plant_names,
                   taxon_name, accepted_plant_name_id, taxon_rank, taxon_status, family, genus, species),
            join_by("scientificName" == "taxon_name"),
            relationship = "many-to-many") %>%
  distinct()

### Pick out whether there is a matching accepted or synonym name
### If there is only one match, keep; 
### If there is more than one match:
###  - Discard if there is no POWO match
###  - Keep if the name is accepted
###  - Discard if the name is a synonym and there is another match which is accepted
###  - Keep if the name is a synonym but there also isn't another match which is accepted
plant_assessment_names <- plant_assessment_names %>%
  group_by(scientificName) %>%
  mutate(powo_id = accepted_plant_name_id) %>%
  arrange(scientificName) %>%
  mutate(keep = case_when(length(accepted_plant_name_id) == 1                              ~ 1, # 72,621 rows
                          length(accepted_plant_name_id)  > 1 & is.na(powo_id)             ~ 0, #  2,370 rows
                          length(accepted_plant_name_id)  > 1 & taxon_status == "Accepted" ~ 1, #  1,977 rows
                          length(accepted_plant_name_id)  > 1 & taxon_status == "Synonym" &  any(taxon_status == "Accepted")  ~ 0, # 168 rows
                          length(accepted_plant_name_id)  > 1 & taxon_status == "Synonym" & !any(taxon_status == "Accepted")  ~ 1, # 147 rows
                          .default = NA
  ))

### Finally, get the final accepted scientific name
plant_assessment_names <- left_join(plant_assessment_names,
                                    select(plant_names, plant_name_id, taxon_name),
                                    join_by("powo_id" == "plant_name_id")) %>%
  rename(powo_name = taxon_name)

taxon_ref_df <- plant_assessment_names %>%
  filter(keep == 1) %>%
  select(scientificName, powo_name, powo_id, taxon_rank, taxon_status) %>%
  mutate(scientificName = gsub(" ", "_", scientificName)) %>%
  mutate(powo_name      = gsub(" ", "_", powo_name))

### 74,735 unique IUCN names
length(unique(taxon_ref_df$scientificName))

### 72,499 unique POWO names
length(unique(taxon_ref_df$powo_name))

### 1,404 IUCN names with no match to POWO
taxon_ref_df %>% filter(is.na(powo_name))

### 843 IUCN names match to multiple POWO names (1,570 individual POWO names)
doubles <- taxon_ref_df$powo_name[duplicated(taxon_ref_df$powo_name)]
taxon_ref_df %>%
  filter(!is.na(powo_name)) %>%
  filter(powo_name %in% doubles[!is.na(doubles)]) %>%
  arrange(powo_name)

### Save to disk
write_csv(taxon_ref_df,
          file.path(threatDir, "plants_iucn_powo_names_table.csv"))
