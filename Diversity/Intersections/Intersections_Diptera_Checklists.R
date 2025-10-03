####################################################################################################
### Run intersections between Diptera (Systema Dipterorum) checklists and bioregions
### Charlie Marsh
### charliem2003@github
### 11/2024
###
### Regions map used: Bioregions_checklists
### Taxonomy used:    Systema Dipterorum
### Checklist used:   Systema Dipterorum
###
### CITATION: Evenhuis, N.L. & Pape, T. (editors). 2024. Systema Dipterorum, Version 5.3.
###   http://diptera.org/ accessed on 27/08/2024
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
library(ggplot2)
library(sf)

### Locations of data, scripts and results
baseDir   <- "/mnt/Work"
regionDir <- file.path(baseDir, "NUS", "BTAS_data", "Bioregions") # dir that contains the regions data
dataDir   <- file.path(baseDir, "spatial_data", "biodiversity", "checklists", "diptera") # dir that contains the checklist data
resDir    <- file.path(baseDir, "NUS", "BTAS", "Intersections") # dir to save results to
funDir    <- file.path(baseDir, "NUS", "BTAS", "Analysis_functions")   # dir that contains the function scripts

### Function for extracting region info from inside brackets
source(file.path(funDir, "Intersections", "extractParantheses.R"))

#==================================================================================================#
#------------------------------------------- Data prep --------------------------------------------#
#==================================================================================================#

### Read in and process data
dip <- read.csv(file.path(dataDir, "Oriental-valid-spp_DIPTERA.csv")) %>%
  as_tibble() %>% 
  mutate(Species = gsub(" ", "_", Valid.Name.Short)) %>%
  filter(Species != "") %>%
  rename(Genus = Valid_Genus_Check) %>%
  select(Species, Genus, Family, Year, Range)

### Global data to get species numbers - 169,674 species
global <- readxl::read_xlsx(file.path(dataDir, "SD-world-valid-extant-species.xlsx"))
length(unique(global$`Valid Name`))

### Spreadsheet with codes to assign to countries and provinces to bioregions
bioregions <- read.csv(file.path(dataDir, "country_codes_BTAS_diptera.csv")) %>%
  rename(Region = Description)

### Spreadsheet outlining bioregions when distribution information is 'X to Y'
paths <- read.csv(file.path(dataDir, "country_paths_BTAS_diptera.csv")) %>%
  select(Raw, Region, AsiaNative, Also_nontrop) %>%
  filter(!is.na(Region) | AsiaNative == 1) %>%
  separate_longer_delim(Region, delim = ", ") %>%
  rename(Bioregion = Region) %>%
  rename(Region = Raw) %>%
  mutate(Exclude = 0)

### Join together to get full list of names
bioregions <- bind_rows(bioregions, paths)

#==================================================================================================#
#------------------------------------ Determine intersections -------------------------------------#
#==================================================================================================#

### Replace square brackets with normal brackets, and semicolons with commas and remove full stops
dip$Range <- gsub("[", "(", dip$Range, fixed = TRUE)
dip$Range <- gsub("]", ")", dip$Range, fixed = TRUE)
dip$Range <- gsub(";", ",", dip$Range, fixed = TRUE)
dip$Range <- gsub(".",  "", dip$Range, fixed = TRUE)

### Make everything lower case for standardisation
dip$Range <- tolower(dip$Range)

### Extract province/island info from brackets and append to country name
dip_region <- dip %>%
  mutate(Range = sapply(Range, extractParantheses))

### Separate out individual entries to unique rows
dip_region <- dip_region %>%
  separate_longer_delim(Range, delim = ", ") %>%
  mutate(Range = gsub(", &", " &", Range)) %>%
  separate_longer_delim(Range, delim = "& ")

### Remove entries for introductions
dip_region <- dip_region %>%
  filter(!grepl("Intro", Range, ignore.case = TRUE, fixed = FALSE))

# ### Remove questionables
# dip <- dip %>%
#   filter(!grepl("[?]", Range, ignore.case = TRUE, fixed = FALSE))

### Remove leading/trailing spaces
dip_region <- dip_region %>%
  mutate(Region = trimws(Range, which = "both"))

### Out of 6420 species with records in China and/or Japan there are 1391 species in raw data with
###   -> retain China and Japan
CN_JP_provinces <- c("indochina", "inochina", "indo-china",
                     "amami", "bonin", "funaura", "hokkaido", "honshu", "iriomote", "ishigaki",
                     "kyushu","okinawa", "omoto-dake", "ryukus", "ryuku", "ryukuys", "ryukyu", 
                     "shikoku", "sikoku", "tsushima", "yaeyama", "yonaha",
                     "anhwei", "anhui", "beijing", "chekiang", "checkiang", "chinghai", "chongqing",
                     "fukian", "fukien", "fujian", "fuijan", "fuijian", "gansu", "guangxi", "gaunxi",
                     "guanxi", "guannxi", "guandong", "guangdong", "guangxa", "guangzhou", "guizhou",
                     "guizhlou", "haidan", "hainan", "hebei", "henan", "honam", "hongkong", "hong kong",
                     "hopeh", "hubei", "hunan", "hupeh", "hunan", "jangxi", "jiang", "jiangxi", "jilin",
                     "kiangsi", "kiangsu", "kwangsi", "kwangting", "kwantung", "kwangtung", "kuling",
                     "liaoning", "macao", "macau", "manchuria", "manshan", "monggol", "neimenggol", "qinghai",
                     "shandong", "shanghai", "shantung", "shaanxi", "shannxi", "shanxi", "sichuan",
                     "sicuan", "sikoku", "szechuan", "szechwan", "china taiwan", "tibet", "tsinghai", "xizang",
                     "xiamen", "yunnan", "yunan", "yunnnan", "yunnsn", "zheijang", "zhejiang", "zuangxi")
length(unique(dip_region$Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), dip_region$Region)]))
dip_region %>%
  filter(Species %in% Species[grepl(paste0("china|japan|", CN_JP_provinces, collapse = "|"), Region)]) %>%
  filter(!Species %in% Species[grepl(paste(CN_JP_provinces, collapse = "|"), Region)]) %>%
  summarise(length(unique(Species)))

### If country == Indonesia or Malaysia (774 species), but more detailed province info also given, remove those rows
MalProvinces <- c("peninsular", "pen", "peninsula", "penisular", "malaysia w", "w ma laysia", "w malaysia",
                  "e malaysia", "east malaysia", "anamba", "anambas", "baha parat","balai ringgin", "bario",
                  "batu niah", "belalong", "bentong", "borneo", "brinchang", "bt timah", "bukit fraser",
                  "bukit frazer", "bukit timah", "cameron highlands", "danum", "genting highlands", "gombak valley",
                  "gunung beremban", "gunung brinchang", "gunung jasar", "hulu langat", "jesselton", "johor",
                  "johore", "k tahan", "kalabakan", "kedah", "kinabalu", "kuala lumpur", "labuan", "malacca",
                  "malatsia", "mendolong", "mengalum", "miri", "mulu", "natuna islands", "negara", "negeri sembilan",
                  "pahang", "papar", "penang", "perak", "rhio", "riau", "ringlet", "sabah", "sandakan",
                  "sarawak", "\nsarawak", "selangor", "selenagor", "semongoh", "simunjan", "singapore",
                  "sipitang", "sungai tekala", "sungei menjuau", "taman negara", "tanah rata", "tenompok",
                  "tioman", "ulu gombak", "ulu langat", "ulu langet")
spMalaysia <- unique(dip_region$Species[dip_region$Species %in% dip_region$Species[dip_region$Region %in% c("malay", "malaysi", "malaysi a", "malaysia","mayalsia", "mayaysia", "mlaya")]])
for(i in 1:length(spMalaysia)) {
  sp <- filter(dip_region, Species == spMalaysia[i])
  if(sum(grepl(paste(MalProvinces,  collapse = "|"), sp$Region)) >= 1) {
    dip_region <- filter(dip_region, Species != spMalaysia[i] | !Region %in% c("malay", "malaysi", "malaysi a", "malaysia","mayalsia", "mayaysia", "mlaya"))
  }
}

IndoProvinces <- c("alor", "amboina", "ambon", "ardjoeno", "aru", "babar",  "bali", "bangka", "batjan",
                   "blume", "bonerate", "borneo", "bukit lawang", "buru", "celebes", "celelbes", "ceram",
                   "dammar", "djakarta", "dumoga-bone", "dyaul", "flores", "gaenaeng gadeh", "gedangan",
                   "gunung leuser", "halmahera", "halmaheira", "iles de la sonde", "irian", "irian jaya",
                   "jambi", "java", "ja va", "ja\\va", "krakatau", "lesser sunda", "lombok", "kai", "kalawaranaputi",
                   "kalimantan", "kangean", "karimundjawa", "kei", "key islands", "kiriwina", "komodo",
                   "malino", "maluku", "madura", "manus", "mataloko", "menado","mentawai", "molucas",
                   "molucca islands", "morotai", "mysol", "new guinea", "new guinra", "new guniea", "new guine a",
                   "new ireland", "nias", "nokilalaki", "nongkodjadjar", "noongan", "nusa ten ggara",
                   "nusa teng gara", "nusa tenggar a", "nusa tenggara",  "papua", "papuan", "palopo",
                   "palu", "pismarck", "png", "pontianak", "sangi", "sangihe", "santong", "sebesi", "sembaloen",
                   "siberut", "simalur", "simuelue", "simueulue", "sipora", "solea", "suban", "sula",
                   "sulawesi", "sulawezi", "dulawesi", "sulaw esi", "sulaw", "sumatera", "sumatra", "sumatria",
                   "sumtra", "sumutra", "sumata", "sumat ra", "soembawa", "sumba", "sumbawa", "tandi andalas",
                   "ternate", "timor", "tondano", "topidi", "tosari", "tjibodas", "wetar", "wonosobo", "uncarya")
spIndonesia <- unique(dip_region$Species[dip_region$Species %in% dip_region$Species[dip_region$Region %in% c("indonesia", "i ndonesia", "indoensia", "indoesia", "indon esia", "indonesai", "oriental region: indonesia")]])
for(i in 1:length(spIndonesia)) {
  sp <- filter(dip_region, Species == spIndonesia[i])
  if(sum(grepl(paste(IndoProvinces,  collapse = "|"), sp$Region)) >= 1) {
    dip_region <- filter(dip_region, Species != spMalaysia[i] | !Region %in% c("indonesia", "i ndonesia", "indoensia", "indoesia", "indon esia", "indonesai", "oriental region: indonesia"))
  }
}

### assign to BTAS bioregions
dip_region <- dip_region %>%
  left_join(bioregions, by = "Region", relationship = "many-to-many") %>%
  filter(Exclude == 0) %>%
  select(Species, Family, Genus, Year, Region, Bioregion, AsiaNative, Also_nontrop)

### to check which countries aren't matched to Bioregions
# sort(unique(select(filter(dip_region, is.na(Bioregion)), Region)$Region))

## which species are native to tropical asia
native <- dip_region %>%
  group_by(Species) %>%
  reframe(AsiaNative   = max(AsiaNative, na.rm = TRUE),
          Also_nontrop = max(Also_nontrop, na.rm = TRUE)) %>%
  mutate(AsiaNative   = case_when(AsiaNative   >= 1 ~ 1,
                                  AsiaNative   == 0 ~ 0)) %>%
  mutate(Also_nontrop = case_when(Also_nontrop >= 1 ~ 1,
                                  Also_nontrop == 0 ~ 0))
dip_region <- dip_region %>%
  select(-AsiaNative, -Also_nontrop, -Region) %>%
  left_join(native, by = "Species") %>%
  # filter(AsiaNative == 1) %>%
  distinct()# %>%
  # mutate(Bioregion = case_when(is.na(Bioregion) & AsiaNative == 1 ~ "Unknown",
  #                              .default = Bioregion))

### which species are endemic to tropical asia
endemic <- dip_region %>%
  group_by(Species) %>%
  reframe(NonEndemic = any(is.na(Bioregion)) | any(Also_nontrop == 1)) %>%
  mutate(AsiaEndemic = case_match(NonEndemic,
                                  FALSE ~ 1,
                                  TRUE  ~ 0))
dip_region <- left_join(dip_region, endemic, by = "Species")

### Remove non-native species and records not assigned to bioregions
dip_region <- dip_region %>%
  filter(AsiaNative == 1) %>%
  filter(!is.na(Bioregion)) %>%
  select(-Also_nontrop, -NonEndemic)

### convert to species x bioregion data frame
dip_region <- dip_region %>%
  filter(!duplicated(.))  %>%
  pivot_wider(id_cols = c(Species, Genus, Family, Year, AsiaEndemic),    
              names_from  = Bioregion,
              values_from = AsiaNative, 
              values_fill = 0) %>%
  mutate(AsiaNative = 1) %>%
  relocate(Species, Genus, Family, Year, AsiaNative, AsiaEndemic)

### Sometimes subspecies are given a unique speciesId so a species may have multiple rows. Merge them
bioregions <- unique(bioregions$Bioregion)[!is.na(unique(bioregions$Bioregion))]
dip_region <- dip_region %>%
  reframe(Genus       = unique(Genus),
          Family      = unique(Family),
          Year        = min(Year),
          AsiaEndemic = min(AsiaEndemic),
          across(any_of(bioregions), ~ max(.)),
          .by = Species)

### We have 300 species where tropical asia distribution is unknown
sum(dip_region$Unknown)

### But there are 112 cases where bioregion is Unknown but the species also has other bioregion info
### In these cases we will remove the Unknown definition
dip_region$Unknown[dip_region$Unknown == 1 & rowSums(dip_region[, bioregions[bioregions != "Unknown"]])] <- 0

### total number of species native to region and number of species endemic to tropical asia
c(native = nrow(dip_region),
  endemic = sum(dip_region$AsiaEndemic))
colSums(dip_region[, bioregions])

### save results
write.csv(dip_region, file.path(resDir, "Intersections", "Intersections_bioregions_diptera.csv"),
          quote = FALSE, row.names = FALSE)

####################################################################################################
### plot richness and turnover

regions <- st_read(file.path(regionDir, "Bioregions_checklists.gpkg")) %>%
  st_simplify(dTolerance = 1000)

rich <- data.frame(Bioregion = names(dip_region[, -c(1:5)]),
                   Rich = colSums(dip_region[, -c(1:5)]))
rich <- left_join(regions, rich, by = "Bioregion")

ggplot() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA),
        # legend.position = c(0.1, 0.02),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.justification = c("left", "bottom"),
        legend.key = element_blank()) +
  scale_fill_viridis_c(option = "C", name = NULL, trans = "log10") +
  guides(size = guide_legend(nrow = 2)) +
  scale_radius(range = c(3, 8),
               breaks = round(10 ^ seq(log10(min(rich$Rich)), log10(max(rich$Rich)), length.out = 5), 0),
               name = "Species Richness",
               trans = "log10") +
  geom_sf(data = rich, aes(fill = Rich)) + 
  geom_sf(data = st_centroid(rich),
          col = "black", pch = 21, aes(size = Rich, fill = Rich))
