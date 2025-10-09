####################################################################################################
### Prepare final bioregions intersections for the threats analysis.
### Charlie Marsh
### charliem2003@github
### 05/2025
###
### This simply takes the intersections for each taxon and:
### 1. Merges Andamans into IndoChina
### 2. Merges Palawan into Philippines
### 3. Removes Northern Australia and species only found there within region
### 4. Removes Northwest and Northeast Australian Shelves and species only found within those regions
###
### The final versions are saved in the parent directory 'Threats'
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Libraries
library(dplyr)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir  <- "Threats"                                             # project dir
divDir   <- file.path("Diversity_and_endemism", "directory")      # path to 'Diversity_and_endemism' folder (for non-IUCN intersections)


projDir <- "/mnt/Work/NUS/Biodiversity-of-tropical-Asia/Threats"
divDir  <- "/mnt/Work/NUS/BTAS"

### You shouldn't need to adjust these folders
iucnDir  <- file.path(projDir, "Intersections", "Intersections", "Final")  # dir with IUCN range map intersections
interDir <- file.path(divDir,  "Intersections", "Intersections")           # dir with intersections in the 'Diversity_and_endemism' folder
resDir   <- file.path(projDir, "Intersections")                            # dir to save results to

### List of taxa
taxa <- read_csv(file.path(projDir, "Global_totals.csv"))
taxa <- taxa$Name

#==================================================================================================#
#------------------------------------------- Processing -------------------------------------------#
#==================================================================================================#

for(taxon in taxa) {
  ### For groups with IUCN range maps get from intersections carried out for threats analysis
  if(taxon %in% c("amphibians", "birds", "freshwater_crabs", "freshwater_fish",  "mammals",
                  "reptiles", "bony_fish", "corals", "sharks")) {
    dd <- read.csv(file.path(iucnDir, paste0("Intersections_bioregions_", taxon, "_5.0_25.csv")))
  }
  
  ### For groups with range maps (but not IUCN) copy from diversity analyses
  if(taxon %in% c("ants")) {
    dd <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", taxon, "_5.0_25.csv")))
  }
  
  ### For other groups with just checklists copy from intersections folder
  if(taxon %in% c("angiosperms", "bees", "butterflies", "chilopods", "diplopods", "diptera", "ferns",
                  "gymnosperms", "lycophytes", "miridae", "phasmids", "spiders", "sponges", "trichoptera")) {
    dd <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", taxon, ".csv")))
  }
  
  ### Merge Andamans into IndoChina
  if(any(names(dd) == "Andamans")) {
    dd <- dd %>%
    mutate(IndoChina = case_when(Andamans == 1 ~ 1,
                                 .default = IndoChina)) %>%
      select(-Andamans)
  }
  
  ### Merge Palawan into Philippines
  if(any(names(dd) == "Palawan")) {
    dd <- dd %>%
      mutate(Philippines = case_when(Palawan == 1 ~ 1,
                                     .default = Philippines)) %>%
      select(-Palawan)
  }
  
  ### Remove australia-only species
  dd <- dd %>%
    select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
    mutate(Sum = rowSums(select(., !any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic", "OutsideAsia"))))) %>%
    filter(Sum > 0) %>%
    select(-Sum)
  
  ### Save to threats directory
  write.csv(dd, file.path(resDir, paste0("Intersections_bioregions_", taxon, ".csv")),
            quote = FALSE, row.names = FALSE)
}

