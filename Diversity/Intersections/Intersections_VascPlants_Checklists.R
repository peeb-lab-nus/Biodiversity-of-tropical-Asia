####################################################################################################
### Calculate checklists for plants and BTAS bioregions using the World Checklist of Vascular
### Plants (WCVP) botanical country checklists. Checklists for provinces of southern China (Hainan, 
### Hong Kong, Guangdong, Guangxi, Macao, Yunnan and Taiwan) have been supplemented with data from
### Xing Yaowu and Harald Schneider.
###
### Lim Jun Ying Charlie Marsh
### charliem2003@github
### 12/2024
###
### Taxonomy used:    WCVP 2024 - downloaded 21/05/2024
### Checklist used:   WCVP 2024 - downloaded 21/05/2024
###
### CITATION:
### Govaerts R (ed.). 2024. WCVP: World Checklist of Vascular Plants. Facilitated by the Royal
###  Botanic Gardens, Kew. [WWW document] URL http://sftp.kew.org/pub/data-repositories/WCVP/
###  [accessed 21 May 2024].
###
### saves csv with the presence-absence of each species within each bioregion. Separate checklists
### are generated for subsets of angiosperms, gymnosperms, lycophytes and ferns.
###
### Notes- Deviations from the other bioregions:
### - Palawan could not individually classified since it is in the PHI BC.
### - Java does not include Bali. Bali is included in Lesser Sunda islands instead
### - Indochina includes species labelled as 'tropical' within relevant parts of mainland China,
###   Taiwan, Hainan or Japan i.e., CHC, CHS, TAI, CHH)
### - Checklists at the northern border would likely include many temperate elements, e.g. Himalayas
### - Northern Australia not included since botanical countries are so large
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

## GENERATE PLANT CHECKLISTS
rm(list = ls())

# PACKAGES ============
library(stringr)
library(vegan)
library(reshape2)
library(plyr)
library(dplyr)
library(readxl)
library(kewr)
library(stringdist)

# DIRECTORIES ========================
if(Sys.info()["user"] == "junyinglim"){
  main.dir <- "~/OneDrive - National University of Singapore/PEEBL Projects/BTAS/"
  powo.dir <- "~/OneDrive - National University of Singapore/PEEBL Research Assets/biodiversity-data/POWO/"
}
if(Sys.info()["user"] == "charlie"){
  base.dir <- file.path("/mnt", "Work")
  main.dir <- file.path(base.dir, "NUS", "BTAS")
  powo.dir <- file.path(base.dir, "spatial_data", "biodiversity", "checklists", "plants", "POWO")
}
misc.dir <- file.path(main.dir, "BTAS_data","taxonomies")
output.dir <- file.path(main.dir, "Intersections", "Intersections")  # dir to save checklists to
fun.dir     <- file.path(main.dir, "Analysis_functions")              # dir that contains the function scripts

### Function for extracting description year from authorships
source(file.path(fun.dir, "Discovery_rates", "clean_publication_dates.R"))

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

# IMPORT DATA ========================
plant_names <- read.table(file.path(powo.dir, "wcvp_v13", "wcvp_names.csv"),
                          header = TRUE, sep = "|", quote = "", fill = TRUE, encoding = "UTF-8")
length(unique(plant_names$plant_name_id))  # 1,431,677
sum(duplicated(plant_names$plant_name_id)) # checking for duplicates

plant_dist <- read.table(file.path(powo.dir, "wcvp_v13","wcvp_distribution.csv"),
                         header = TRUE, sep = "|", quote = "", fill = TRUE, encoding = "UTF-8")
length(unique(plant_dist$plant_name_id))   # 440,341 names
all(unique(plant_dist$plant_name_id) %in% unique(plant_names$plant_name_id)) # all are in plant_name

# Are all names in the distribution, also accepted names? No, but mostly
head(plant_names)
table(subset(plant_names, plant_name_id %in% unique(plant_dist$plant_name_id))$taxon_status)

#==================================================================================================#
#--------------------------------- Extract tropical China species ---------------------------------#
#==================================================================================================#

## OF THE SPECIES IN THE SOUTH OF CHINA, HOW COMPLETE IS THE CLIMATE ASSOCIATION DATA ===========
tdwg_trop_asia_other <- c("ASS", "BAN", "EHM", "IND", "NEP", "SRL", "WHM", # Indian subcontinent
                          "LDV", "MDV", # Laccadives and Maldives,
                          "AND", "NCB", # Andaman and Nicobar
                          "CBD", "LAO", "MYA", "THA", "VIE", "SCS", # Indochina
                          "BOR", # Borneo
                          "JAW", # Jawa
                          "LSI", # Lesser sundas
                          "MLY", # Malaya
                          "MOL", # Maluku
                          "PHI", # Philippines
                          "SUL", # Sulawesi
                          "SUM", # Sumatra
                          "BIS", "NWG") # New Guinea
# "XMS", "CKI",  "SOL"
# (Christmas isl., Cocos-Keeling Isl., Solomons excluded)
# "PAK" (Pakistan)

trop_asia_name_ids <- plant_dist %>%
  subset(introduced == 0 & location_doubtful == 0) %>%
  subset(area_code_l3 %in% tdwg_trop_asia_other) %>%
  left_join(plant_names, by = "plant_name_id") %>%
  subset(taxon_rank == "Species" & taxon_status == "Accepted") %>%
  .[["plant_name_id"]] %>%
  unique() # list of names excluding china

china_name_ids <- plant_dist %>% 
  subset(introduced == 0 & location_doubtful == 0) %>%
  subset(area_code_l3 %in% c("CHC", "CHS", "TAI", "CHH", "NNS")) %>%
  left_join(plant_names, by = "plant_name_id") %>%
  subset(taxon_rank == "Species" & taxon_status == "Accepted") %>%
  .[["plant_name_id"]] %>%
  unique()
china_names <- subset(plant_names, plant_name_id %in% china_name_ids)

length(china_names$plant_name_id[grepl(china_names$climate_description, pattern = "tropical")]) # 13201

table(china_names$climate_description)
temp <- subset(china_names, climate_description == "")
table(temp$climate_description)
sum(temp$plant_name_id %in% trop_asia_name_ids) # 188 are in the rest of tropical asia (which we should include)
# Note: There are non-tropical elements that in other parts of tropical asia,
# so you will get more than 188. From a consistency point of view, they have to be included.
subset(china_names, (!grepl(climate_description, pattern = "tropical")) & 
         plant_name_id %in% trop_asia_name_ids) %>%
  nrow() # 3605
# 3605 + 13201 = 16806

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

# CREATE GLOBAL CHECKLIST ========================
# Clean up plant distribution data
# Only included species rank and accepted taxon status
plant_dist_clean <- plant_dist %>%
  subset(introduced == 0 & location_doubtful == 0) %>%
  select(plant_name_id, area_code_l3) %>%
  left_join(plant_names[c("plant_name_id",
                          "taxon_rank",
                          "taxon_status",
                          "basionym_plant_name_id",
                          "first_published",
                          "family",
                          "genus",
                          "species",
                          "climate_description")],
            by = "plant_name_id") %>%
  mutate(scientific_name = paste(genus, species, sep = " ")) %>%
  subset(taxon_rank == "Species" & taxon_status == "Accepted")

# Assign botanical countries to BTAS bioregions
plant_dist_clean$region <- NA
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("WHM", "NEP",
                                                             "IND", "EHM",
                                                             "ASS", "BAN",
                                                             "SRL")]        <- "Indian_Subcontinent"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("LDV", "MDV")] <- "Indian_Subcontinent"# "Laccadives and Maldives"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("MYA", "LAO",
                                                             "VIE", "THA",
                                                             "CBD", "SCS")] <- "IndoChina"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("CHC", "CHS", "NNS", "CHH", "TAI") &
                        (grepl(plant_dist_clean$climate_description, pattern = "tropical") | 
                               plant_dist_clean$plant_name_id %in% trop_asia_name_ids)] <- "IndoChina"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("MLY")]        <- "Malaya"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("SUM")]        <- "Sumatra"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("PHI")]        <- "Philippines"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("JAW")]        <- "Java"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("SUL")]        <- "Sulawesi"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("BOR")]        <- "Borneo"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("LSI")]        <- "Lesser_Sundas"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("AND", "NCB")] <- "Andamans"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("MOL")]        <- "Maluku"
plant_dist_clean$region[plant_dist_clean$area_code_l3 %in% c("NWG", "BIS")] <- "New_Guinea"

# List of botanical countries included in BTAS bioregions
trop_asia_bcs <- unique(plant_dist_clean$area_code_l3[!is.na(plant_dist_clean$region)])

## CALCULATE ENDEMISM AND TURNOVER ==============================
# 74,662 species in "tropical asia"
length(unique(plant_dist_clean$plant_name_id[!is.na(plant_dist_clean$region)])) 

plant_endemism <- ddply(.data = plant_dist_clean,
                        .variables = .(plant_name_id),
                        .fun = summarise,
                        AsiaEndemic = all(area_code_l3 %in% trop_asia_bcs),
                        .progress = "text")

sum(plant_endemism$AsiaEndemic) # 70421 species are endemic to "tropical asia"

## GENERATE REGIONAL CHECKLIST ==============================
trop_asia_plant_mat <- acast(formula = plant_name_id ~ region,
                             fill = 0, value.var = "region",
                             fun.aggregate = function(x){length(unique)},
                             data = subset(plant_dist_clean, !is.na(region))) %>%
  as.data.frame() %>%
  mutate(plant_name_id = as.numeric(rownames(.))) %>%
  left_join(plant_names[c("plant_name_id",
                          "taxon_rank",
                          "taxon_status",
                          "basionym_plant_name_id",
                          "first_published",
                          "family",
                          "genus",
                          "species",
                          "taxon_name")],
            by = "plant_name_id")

# clean_publication_dates <- function(x, max = TRUE){
#   # Cleaning publication dates from the WCVP data
#   #
#   # Arguments:
#   #   x, chr, vector of publication dates (`first_published`)  
#   #   max, logical, whether to choose the oldest year in the string (if FALSE, then minimum is taken)
#   #
#   # Returns:
#   #   numeric
#   
#   res <- vector()
#   
#   for(i in 1:length(x)){
#     if(x[i] == ""){
#       res[i] <- NA
#     } else {
#       match_positions <- gregexpr("[0-9]+", text = x[i])
#       numbers <- regmatches(x[i], match_positions)  
#       numbers_vector <- as.numeric(unlist(numbers))
#       if(max == TRUE){
#         res[i] <- max(numbers_vector)
#       } else {
#         res[i] <- min(numbers_vector)
#       }  
#     }
#   }
#   return(res)
# }

nrow(trop_asia_plant_mat) # 74,662
sum(trop_asia_plant_mat$first_published == "") # 45 missing data points
sum(is.na(trop_asia_plant_mat$first_published)) # no NAs
sum(is.na(trop_asia_plant_mat$basionym_plant_name_id)) # 53211 out of 67186

# Extract publication years
trop_asia_plant_mat2 <- trop_asia_plant_mat %>% 
  dplyr::rename("taxon_published_year" = "first_published") %>%
  mutate("taxon_published_year" = clean_publication_dates(taxon_published_year))

# Not all taxa have publication year - try and get from basionym instead
trop_asia_plant_mat2$basionym_ref <- ifelse(is.na(trop_asia_plant_mat2$basionym_plant_name_id),
                                            trop_asia_plant_mat2$plant_name_id,
                                            trop_asia_plant_mat2$basionym_plant_name_id)

# Join everything together
trop_asia_plant_mat_final <- trop_asia_plant_mat2 %>%
  left_join(plant_endemism[c("plant_name_id", "AsiaEndemic")]) %>%
  mutate(AsiaEndemic = as.numeric(AsiaEndemic)) %>%
  left_join(plant_names[c("plant_name_id", "first_published")],
            by = c("basionym_ref" = "plant_name_id")) %>%
  mutate("basionym_published_year" = clean_publication_dates(first_published)) %>%
  rename(Species = species) %>%
  rename(Genus = genus) %>%
  mutate(Species = paste(Genus, Species, sep = "_")) %>%
  mutate(Family = family) %>%
  select(c("plant_name_id",
           "taxon_rank",
           "taxon_status",
           "basionym_plant_name_id",
           "taxon_published_year",
           "basionym_published_year",
           "AsiaEndemic",
           "Family",
           "Genus",
           "Species",
           "taxon_name",
           "Indian_Subcontinent",
           # "Laccadives and Maldives",
           "IndoChina",
           "Borneo",
           "Malaya",
           "Philippines",
           "Sumatra",
           "Java",
           "Maluku",
           "Sulawesi",
           "Lesser_Sundas",
           "New_Guinea",
           "Andamans"))

sum(is.na(trop_asia_plant_mat_final$taxon_published_year))          # 45 species without publication year
sum(is.na(trop_asia_plant_mat_final$basionym_published_year))       # 2 species without basionym year
range(trop_asia_plant_mat_final$taxon_published_year, na.rm = TRUE) # 1753 to 2024

# Final data frame
trop_asia_plant_mat_final <- trop_asia_plant_mat_final %>%
  rename(Year = basionym_published_year) %>%
  select(Species, Genus, Family, Year, AsiaEndemic,
         Indian_Subcontinent, IndoChina, Andamans, Philippines, Malaya, Sumatra, Java, Borneo,
         Sulawesi, Lesser_Sundas, Maluku, New_Guinea)

#==================================================================================================#
#---------------------------------- Subset plant groups and save ----------------------------------#
#==================================================================================================#

all_plant_families <- unique(trop_asia_plant_mat_final$Family)

# Ferns
fern_families <- c("Gleicheniaceae",
                   "Cyatheaceae",
                   "Polypodiaceae",
                   "Schizaeaceae",
                   "Osmundaceae",
                   "Equisetaceae",
                   "Pteridaceae",
                   "Aspleniaceae",
                   "Marattiaceae",
                   "Ophioglossaceae",
                   "Psilotaceae",
                   "Hymenophyllaceae",
                   "Dipteridaceae",
                   "Matoniaceae",
                   "Marsileaceae",
                   "Dennstaedtiaceae",
                   "Salviniaceae",
                   "Cystodiaceae",
                   "Saccolomataceae")
length(unique(plant_dist_clean$plant_name_id[plant_dist_clean$family %in% fern_families]))      # 11,885 global species                                          # 3,903 species
sum(trop_asia_plant_mat_final$Family %in% fern_families)                                        #  3,903 species
sum(trop_asia_plant_mat_final$AsiaEndemic[trop_asia_plant_mat_final$Family %in% fern_families]) #  2,749 endemic species

# Lycophytes
lycophyte_families <- c("Lycopodiaceae",
                        "Selaginellaceae",
                        "Isoetaceae")
length(unique(plant_dist_clean$plant_name_id[plant_dist_clean$family %in% lycophyte_families]))      # 1,438 global species                                          # 3,903 species
sum(trop_asia_plant_mat_final$Family %in% lycophyte_families)                                        #   389 species
sum(trop_asia_plant_mat_final$AsiaEndemic[trop_asia_plant_mat_final$Family %in% lycophyte_families]) #   302 endemic species

# Gymnosperms
gymnosperm_families <- c("Araucariaceae",
                         "Cupressaceae",
                         "Pinaceae",
                         "Podocarpaceae",
                         "Taxaceae",
                         "Ephedraceae",
                         "Gnetaceae",
                         "Cycadaceae",
                         "Cephalotaxaceae")
length(unique(plant_dist_clean$plant_name_id[plant_dist_clean$family %in% gymnosperm_families]))      # 946 global species                                          # 3,903 species
sum(trop_asia_plant_mat_final$Family %in% gymnosperm_families)                                        # 288 species
sum(trop_asia_plant_mat_final$AsiaEndemic[trop_asia_plant_mat_final$Family %in% gymnosperm_families]) # 266 endemic species

# Angiosperms (everything else)
angiosperm_families <- all_plant_families[!all_plant_families %in% c(gymnosperm_families,
                                                                     lycophyte_families,
                                                                     fern_families)]
length(unique(plant_dist_clean$plant_name_id[plant_dist_clean$family %in% angiosperm_families]))      # 338,538 global species                                          # 3,903 species
sum(trop_asia_plant_mat_final$Family %in% angiosperm_families)                                        #  70,082 species
sum(trop_asia_plant_mat_final$AsiaEndemic[trop_asia_plant_mat_final$Family %in% angiosperm_families]) #  59,631 endemic species

# Save csvs
# write.csv(trop_asia_plant_mat_final,
#           file.path(output.dir, "Intersections_bioregions_vascularplants.csv"),
#           row.names = FALSE)

write.csv(subset(trop_asia_plant_mat_final, Family %in% lycophyte_families),
          file.path(output.dir, "Intersections_bioregions_lycophytes.csv"),
          row.names = FALSE)

write.csv(subset(trop_asia_plant_mat_final, Family %in% fern_families),
          file.path(output.dir, "Intersections_bioregions_ferns.csv"),
          row.names = FALSE)

write.csv(subset(trop_asia_plant_mat_final, Family %in% gymnosperm_families),
          file.path(output.dir, "Intersections_bioregions_gymnosperms.csv"),
          row.names = FALSE)

write.csv(subset(trop_asia_plant_mat_final, Family %in% angiosperm_families),
          file.path(output.dir, "Intersections_bioregions_angiosperms.csv"),
          row.names = FALSE)

