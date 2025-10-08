####################################################################################################
### Intersections between Chilopods and bioregions
### Chilobase 2.0 integrated and updated and corrected based on the scientific literature published
### between 2000 and 2023 by Lucio Bonato in 2024
###
### Charlie Marsh
### charliem2003@github
### 12/2024
###
### Taxonomy used:    https://chilobase.biologia.unipd.it/
###                   GBIF taxonomic backbone is needed for family info (https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c)
### Checklist used:   https://chilobase.biologia.unipd.it/
###
### CITATION:
### Bonato L., Chagas Junior A., Edgecombe G.D. Lewis J.G.E., Minelli A., Pereira L.A., Shelley R.M.,
###  Stoev P., Zapparoli M. (2016) ChiloBase 2.0 - A World Catalogue of Centipedes (Chilopoda).
###  Available at https://chilobase.biologia.unipd.it. 
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
library(stringr)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir   <- "Diversity_and_endemism"                    # project dir
dataDir  <- file.path("ChiloBase", "data", "directory")  # dir that contains the checklist data
taxonDir <- file.path("GBIF", "backbone", "directory")   # dir that contains GBIF taxonomic backbone (for families)

### You shouldn't need to adjust these folders
resDir   <- file.path(projDir, "Intersections")          # dir to save results to

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Read in cleaned ChiloBase data - 488 species (global richness = 3340 in ChiloBase 2.0)
all <- read.csv(file.path(dataDir, "Tropical Asia - Chilopoda - by Lucio Bonato.csv")) %>%
  mutate(species = trimws(species, "both")) %>%
  as_tibble()

### Taxonomic information with authorship details
taxon <- read.csv(file.path(dataDir, "Tropical Asia - Chilopoda - by Lucio Bonato - taxonomy.csv")) %>%
  select(full.species.name, genus.name, year)

### There are some weird space punctutation that needs replacing in order to get name matches
taxon <- taxon %>%
  mutate(full.species.name = gsub(strsplit(full.species.name[25], "")[[1]][8], " ", full.species.name)) %>%
  mutate(full.species.name = trimws(full.species.name, "both")) %>%
  mutate(genus.name = trimws(genus.name, "both"))
         
### The taxon info doesn't include families, pull those from the gbif taxonomic backbone
gbif <- read_delim(file.path(taxonDir, "Taxon.tsv"), delim = "\t") %>%
  select(family, genus) %>% 
  distinct() %>%
  filter(!is.na(family)) %>%
  filter(!is.na(genus))

taxon <- taxon %>%
  left_join(gbif,
            join_by("genus.name" == "genus"),
            relationship = "many-to-many")

### The family of Brachigeophilus robustus is unknown
taxon$family[taxon$full.species.name == '"Brachigeophilus" robustus'] <- NA

### There are also some duplicated species with multiple species names
taxon <- taxon %>%
  mutate(family = case_when(family == "Brachyceridae" ~ "Cryptopidae", # There is also a weevil genus Crytops in family Brachyceridae
                            family == "Columbidae"    ~ "Geophilidae", # There is a pigeon genus Geophilus in family Columbidae
                            .default = family)) %>%
  distinct()

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### Pivot columns to longer format
all <- all %>%
  pivot_longer(cols = !"species", names_to = "Region", values_to = "Presence") %>%
  filter(!is.na(Presence))

### Merge with taxonomic information
all <- all %>%
  left_join(taxon, join_by("species" == "full.species.name")) %>%
  mutate(Species = gsub(" ", "_", species)) %>%
  dplyr::rename(Genus   = genus.name) %>%
  dplyr::rename(Year    = year) %>%
  dplyr::rename(Family  = family)

### There is one species where we need to add details manually based on https://eol.org/pages/55690986
all$Family[all$species == "Orphnaeus aporus"] <- "Oryidae"
all$Genus [all$species == "Orphnaeus aporus"] <- "Orphnaeus"
all$Year  [all$species == "Orphnaeus aporus"] <- 1930

### Assign ChiloBase regions to BTAS bioregions
all <- all %>%
  mutate(Bioregion = case_when(Region %in% c("Nepal..NP.", "Bhutan..BT.", "India..IN.",
                                             "Sri.Lanka..LK.", "Maldives..MV.",
                                             "Laccadive.Is...IN.", "Bangladesh..BD.")         ~ "Indian_Subcontinent",
                               Region %in% c("Andaman.Is...IN.", "Nicobar.Is...IN.")          ~ "Andamans",
                               Region %in% c("Myanmar..MM.", "Yunnan", "Guangxi.Zhuang",
                                             "Guangdong", "Hong.Kong..HK.", "Hainan",
                                             "Thailand..TH.", "Laos..LA.", "Vietnam..VN.",
                                             "Cambodia..KH.", "Southern.Ryukyu.Ids..JP.",
                                             "Taiwan..TW.")                                   ~ "IndoChina",
                               Region %in% c("Peninsular.Malaysia..MY.", "Singapore..SG.")    ~ "Malaya",
                               Region %in% c("Philippines..PH.")                              ~ "Philippines",
                               Region %in% c("Sumatera..ID.")                                 ~ "Sumatra",
                               Region %in% c("Jawa..ID.", "Bali..ID.")                        ~ "Java",
                               Region %in% c("Borneo..BOR.")                                  ~ "Borneo",
                               Region %in% c("Sulawesi..ID.")                                 ~ "Sulawesi",
                               Region %in% c("Lesser.Sunda.Islands")                          ~ "Lesser_Sundas",
                               Region %in% c("Maluku..ID.")                                   ~ "Maluku",
                               Region %in% c("Irian.Jaya..ID.", "Papua.New.Guinea..PG.",
                                             "Bismarck.Archipelago..PG.")                     ~ "New_Guinea",
                               Region %in% c("Cocos..Keeling..Is...CC.", "Christmas.I...CX.",
                                             "OUTSIDE.SOUTH.EAST.ASIA")                       ~ "OutsideAsia"))

### Pivot to wide format
all <- all %>%
  select(Species, Genus, Family, Year, Presence, Bioregion) %>%
  filter(Presence == 1) %>%
  distinct() %>%
  pivot_wider(id_cols     = c(Species, Genus, Family, Year),
              names_from  = Bioregion,
              values_from = Presence,
              values_fill = 0)

### Assign endemism
all <- all %>%
  mutate(AsiaEndemic = case_when(OutsideAsia == 0 ~ 1,
                                 OutsideAsia == 1 ~ 0)) %>%
  select(-OutsideAsia) %>%
  relocate(AsiaEndemic, .after = Year)

### Check there are no non-Asian species still remaining
any(rowSums(all[, 6:17]) == 0)

### total number of species native to region and number of species endemic to tropical asia
table(all$AsiaEndemic)

### save results
write.csv(all, file.path(resDir, "Intersections", "Intersections_bioregions_chilopods.csv"),
          quote = FALSE, row.names = FALSE)

#==================================================================================================#
#----------------------------------------- Clean up memory ----------------------------------------#
#==================================================================================================#

rm(list = ls())
