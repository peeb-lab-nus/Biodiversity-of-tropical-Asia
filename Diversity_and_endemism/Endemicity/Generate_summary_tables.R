####################################################################################################
### 
### Charlie Marsh
### charliem2003@github
### 12/2024
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

### Libraries
library(readr)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)

# Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir <- "Diversity_and_endemism"

# You shouldn't need to adjust these folders
interDir <- file.path(projDir, "Intersections", "Intersections")  # dir with intersections data
resDir   <- file.path(projDir, "Endemicity", "Tables")            # dir to save results to
if(!dir.exists(resDir)) { dir.create(resDir, recursive = TRUE) }

### Spreadsheet with taxon information
info <- read.csv(file.path(projDir, "Global_totals.csv")) %>%
  arrange(desc(Group), desc(Taxon_Higher), Taxon)

#==================================================================================================#
#---------------------------------- Summary for all tropical Asia ---------------------------------#
#==================================================================================================#

### List of taxa to plot and bioregion range map type for plotting
taxon_list <- info$Taxon

ddAll <- tibble()
for(taxon in taxon_list) {
  ### Extract taxon info from global info
  ddType    <- info$Type[info$Taxon == taxon]
  ddName    <- info$Name[info$Taxon == taxon]
  totalRich <- info$Global_total[info$Taxon == taxon]
  
  ### For checklists we just need the single csv file
  if(ddType %in% c("Checklist", "Occurrence points")) {
    pa <- read_csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, ".csv")),
                   show_col_types = FALSE, progress = FALSE)
    
    ### Remove australia and any species only found there
    pa <- pa %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic", "OutsideAsia", "AsiaNative")))) > 0)
    
    ddTaxon <- tibble("Taxonomic_group"           = taxon,
                      "No._of_species"            = as.character(format(nrow(pa), big.mark=",")),
                      "Prop._of_global_diversity" = as.character(round(nrow(pa) / totalRich, 3)),
                      "No._of_endemic_species"    = as.character(format(sum(pa$AsiaEndemic), big.mark=",")),
                      "Prop._of_species_endemic"  = as.character(round(sum(pa$AsiaEndemic) / nrow(pa), 3)),
                      "Data_type"                 = info$Type[info$Taxon == taxon],
                      "Data_source"               = info$Data_sources[info$Taxon == taxon])
  }
  
  ### For range maps we need the mean, upper and lower estimates
  if(ddType == "Range maps") {
    paLow  <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, "_7.5_35.csv")))
    paMean <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, "_5.0_25.csv")))
    paUpp  <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, "_2.5_15.csv")))
    
    ### Remove australia and any species only found there
    paLow <- paLow %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic", "OutsideAsia")))) > 0)
    paMean <- paMean %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic", "OutsideAsia")))) > 0)
    paUpp <- paUpp %>%
      select(-any_of(c("Northern_Australia", "Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
      filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic", "OutsideAsia")))) > 0)
    
    ddTaxon <- tibble("Taxonomic_group"           = taxon,
                      "No._of_species"            = paste0(format(nrow(paMean), big.mark=","), "\n", "(",
                                                           format(nrow(paLow), big.mark=","), " - ",
                                                           format(nrow(paUpp), big.mark=","), ")"),
                      "Prop._of_global_diversity" = paste0(round(nrow(paMean) / totalRich, 3), "\n", "(",
                                                           round(nrow(paLow)  / totalRich, 3), " - ",
                                                           round(nrow(paUpp)  / totalRich, 3), ")"),
                      "No._of_endemic_species"    = paste0(format(sum(paMean$AsiaEndemic), big.mark=","), "\n", "(",
                                                           format(sum(paLow$AsiaEndemic), big.mark=","), " - ",
                                                           format(sum(paUpp$AsiaEndemic), big.mark=","), ")"),
                      "Prop._of_species_endemic"  = paste0(round(sum(paMean$AsiaEndemic) / nrow(paMean), 3), "\n", "(",
                                                           round(sum(paLow$AsiaEndemic)  / nrow(paLow), 3), " - ",
                                                           round(sum(paUpp$AsiaEndemic)  / nrow(paUpp), 3), ")"),
                      "Data_type"                 = info$Type[info$Taxon == taxon],
                      "Data_source"               = info$Data_sources[info$Taxon == taxon])
  }
  print(ddTaxon)
  
  ddAll <- bind_rows(ddAll, ddTaxon)
  rm(ddTaxon)
}
write.csv2(ddAll, file.path(resDir, "Trop_Asia_summary.csv"), eol = "\r\n", row.names = FALSE)

#==================================================================================================#
#------------------------------------- Summary for bioregions -------------------------------------#
#==================================================================================================#

### List of taxa to plot and bioregion range map type for plotting
taxon_list <- info$Taxon[info$Group != "Marine"]

### Order of bioregions
bioregions <- c("Indian_Subcontinent", "IndoChina", "Philippines", "Malaya", "Sumatra",
                "Borneo", "Java", "Sulawesi", "Lesser_Sundas", "Maluku", "New_Guinea")

ddAll <- tibble()
for(taxon in taxon_list) {
  ### Extract taxon info from global info
  ddType    <- info$Type[info$Taxon == taxon]
  ddName    <- info$Name[info$Taxon == taxon]
  
  ### For range maps just stick to the mean estimate
  if(ddType %in% c("Checklist", "Occurrence points")) {
    pa <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, ".csv")))
  }
  if(ddType == "Range maps") {
    pa <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, "_5.0_25.csv")))
  }
  
  ### Remove australia and any species only found there
  pa <- pa %>%
    select(-any_of("Northern_Australia")) %>%
    filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic")))) > 0)
  
  ### Merge Palawan and Philippines
  if(any(names(pa) == "Palawan")) {
    pa <- pa %>%
      mutate(Philippines = case_when(Philippines == 0 & Palawan == 0 ~ 0,
                                     Philippines == 1 | Palawan == 1 ~ 1)) %>%
      select(-Palawan)
  }
  
  ### Merge Andamans into IndoChina
  if(any(names(pa) == "Andamans")) {
    pa <- pa %>%
      mutate(IndoChina = case_when(IndoChina == 0 & Andamans == 0 ~ 0,
                                   IndoChina == 1 | Andamans == 1 ~ 1)) %>%
      select(-Andamans)
  }
  
  ### Generate new data frame for if species is bioregion endemic
  endemic <- pa %>%
    mutate(endemic = case_when(rowSums(.[, bioregions]) == 1 & AsiaEndemic == 1 ~ 1,
                               rowSums(.[, bioregions]) == 1 & AsiaEndemic == 0 ~ 0,
                               rowSums(.[, bioregions])  > 1 ~ 0)) %>%
    filter(endemic == 1) %>%
    select(all_of(bioregions))
  
  ### Order bioregions
  pa <- pa %>%
    select(all_of(bioregions))
  
  ### Summarise
  ddTaxon <- c(taxon, paste0(format(colSums(pa), big.mark=","), "\n",
                             "(", round(colSums(endemic) / colSums(pa), 2), ")"))
  names(ddTaxon) <- c("Taxonomic_group", bioregions)
  
  ddAll <- bind_rows(ddAll, ddTaxon)
  rm(ddTaxon)
}

### Tidy up names etc
ddAll <- ddAll %>%
  dplyr::rename(South_Asia           = Indian_Subcontinent,
                Indochina            = IndoChina,
                Malaya_Peninsula     = Malaya,
                Lesser_Sunda_Islands = Lesser_Sundas) %>%
  relocate(c(Taxonomic_group, South_Asia, Indochina, Philippines, Borneo, Malaya_Peninsula, 
             Sumatra, Java, Sulawesi, Lesser_Sunda_Islands, Maluku, New_Guinea))

write.csv2(ddAll, file.path(resDir, "Bioregions_summary.csv"), eol = "\r\n", row.names = FALSE)

#==================================================================================================#
#---------------------------------- Summary for marine ecoregions----------------------------------#
#==================================================================================================#

### List of taxa to plot and bioregion range map type for plotting
taxon_list <- info$Taxon[info$Group == "Marine"]

### Order of bioregions
bioregions <- c("Central_Indian_Ocean_Islands", "West_and_South_Indian_Shelf", "Bay_of_Bengal", "Andaman",
                "Sunda_Shelf", "Java_Transitional", "South_China_Sea", "Western_Coral_Triangle",
                "South_Kuroshio", "Sahul_Shelf", "Eastern_Coral_Triangle")

ddAll <- tibble()
for(taxon in taxon_list) {
  ### Extract taxon info from global info
  ddType    <- info$Type[info$Taxon == taxon]
  ddName    <- info$Name[info$Taxon == taxon]
  
  ### For range maps just stick to the mean estimate
  if(ddType %in% c("Checklist", "Occurrence points")) {
    pa <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, ".csv")))
  }
  if(ddType == "Range maps") {
    pa <- read.csv(file.path(interDir, paste0("Intersections_bioregions_", ddName, "_5.0_25.csv")))
  }
  
  ### Remove australia and any species only found there
  pa <- pa %>%
    select(-any_of(c("Northwest_Australian_Shelf", "Northeast_Australian_Shelf"))) %>%
    filter(rowSums(select(., -any_of(c("Species", "Genus", "Family", "Year", "Range_area", "AsiaEndemic")))) > 0)
  
  ### Generate new data frame for if species is bioregion endemic
  endemic <- pa %>%
    mutate(endemic = case_when(rowSums(.[, bioregions]) == 1 & AsiaEndemic == 1 ~ 1,
                               rowSums(.[, bioregions]) == 1 & AsiaEndemic == 0 ~ 0,
                               rowSums(.[, bioregions])  > 1 ~ 0)) %>%
    filter(endemic == 1) %>%
    select(all_of(bioregions))
  
  ### Order bioregions
  pa <- pa %>%
    select(all_of(bioregions))
  
  ddTaxon <- c(taxon, paste0(format(colSums(pa), big.mark=","), "\n",
                             "(", round(colSums(endemic) / colSums(pa), 2), ")"))
  names(ddTaxon) <- c("Taxonomic_group", bioregions)
  
  ddAll <- bind_rows(ddAll, ddTaxon)
  rm(ddTaxon)
}
write.csv2(ddAll, file.path(resDir, "Marine_ecoregions_summary.csv"), eol = "\r\n", row.names = FALSE)
