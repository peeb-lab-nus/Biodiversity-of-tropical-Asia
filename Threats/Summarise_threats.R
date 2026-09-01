####################################################################################################
### Summarises IUCN threat assessments
### Charlie Marsh
### charliem2003@github
### 07/2026
###
### Human Footprint calculated by Mu et al. (2022) A global record of annual terrestrial Human
### Footprint dataset from 2000 to 2018. Sci Data 9, 176 is available from:
### https://www.nature.com/articles/s41597-022-01284-8.
### Data from 2018 and 2024 needs to be obtained from the authors. Once obtained summarise HFP for
### each year within all cells within each bioregion, and generate a spreadsheet containing the
### columns: 'bioregion', 'year', mean_hfp', 'median_hfp', 'q25_hfp' and 'q75_hfp'. Place in dataDir
###
####################################################################################################

#==================================================================================================#
#--------------------------------------------- Set-up ---------------------------------------------#
#==================================================================================================#

rm(list = ls())

### Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(scales)
library(lubridate)
library(ggplot2)
library(scatterpie)
library(patchwork)
library(sf)

### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir <- "Threats"                                   # project dir

### You shouldn't need to adjust these folders
dataDir <- file.path(projDir, "Data")                  # dir with data and maps
resDir  <- file.path(projDir, "Results")               # dir with summarised assessments
figDir  <- file.path(projDir, "Figures")               # dir to save figures to

#==================================================================================================#
#----------------------------------------- Read in data -------------------------------------------#
#==================================================================================================#

### Taxon group information
taxon_info  <- read_csv(file.path(projDir, "Global_totals.csv")) %>%
  filter(Taxon != "Sponges") %>%
  select(Taxon, Group, Taxon_Higher, Type, Name, Colour) %>%
  mutate(Group = gsub(" NA", "", paste(Group, Taxon_Higher, sep = " "))) %>%
  mutate(GroupForPlot = gsub(" ", "\n", Group, fixed = TRUE)) %>%
  mutate(GroupForPlot = gsub("and\n", "& ", GroupForPlot, fixed = TRUE)) %>%
  mutate(TaxonForPlot = case_when(Name == "phasmids" ~ "Stick and leaf insects",
                                  Name == "sharks"   ~ "Sharks and rays",
                                  .default = Taxon))

### Bioregion information
region_info <- read_csv(file.path(projDir, "Bioregion_info.csv")) %>%
  mutate(LabelForPlot = gsub(" ", "\n", Label, fixed = TRUE))

regions <- region_info %>%
  filter(!is.na(Order)) %>%
  pull(LabelForPlot)

### HFP info from Mu et al 2022 -EDIT FOR YOUR FILE NAME
hfp <- read.csv(file.path(dataDir, "YOUR_HFP_SUMMARY_RESULTS.csv")) %>%
  left_join(region_info, by = "Bioregion") %>%
  mutate(LabelForPlot = factor(LabelForPlot,
                               levels = regions))

### Intersections for all tropical Asia species
all_intersections <- read_csv(file.path(resDir, "all_intersections.csv"))

### All assessments - 176,741 rows for 147,253 species
all_assessments <- read_csv(file.path(resDir, "all_assessments.csv"))
nrow(all_assessments)
length(unique(all_assessments$taxon_scientific_name))

### BTAS map for plotting
base_map        <- readRDS(file.path(dataDir, "gadm_base_map.rds"))
terr_bioregions <- readRDS(file.path(dataDir, "terr_bioregions.rds"))
mare_bioregions <- readRDS(file.path(dataDir, "mare_bioregions.rds")) %>%
  st_crop(st_bbox(base_map))

### Get dividing lines between subregions with land-sea border removed
base_lines <- st_union(st_cast(select(base_map, geom), "MULTILINESTRING"))
mare_lines <- st_union(st_cast(select(mare_bioregions, geom), "MULTILINESTRING"))
terr_lines <- st_union(st_cast(select(terr_bioregions, geom), "MULTILINESTRING"))

base_map_buffer <- st_buffer(base_lines, 10000)
terr_divides <- st_difference(terr_lines, base_map_buffer)
mare_divides <- st_difference(mare_lines, base_map_buffer)

### Colour scheme for Red List categories
colsRedList <- c("CR"             = "#A91F13",
                 "EN"             = "orange2",
                 "VU"             = "#F9DA14",
                 "Non-threatened" = "grey80",
                 "DD"             = "grey60",
                 "Not evaluated"  = "grey40")

####################################################################################################
### Join assessment information to distributions
####################################################################################################

### Tidy up RL cat names and remove species with no assessments
all_assessments <- all_assessments %>%
  mutate(red_list_category_code = case_when(
    red_list_category_code %in% c("CR", "CT")                                  ~ "CR",
    red_list_category_code %in% c("EN")                                        ~ "EN",
    red_list_category_code %in% c("VU", "V")                                   ~ "VU",
    red_list_category_code %in% c("LC", "LR/lc", "NT", "LR/nt", "LR/cd", "NR") ~ "Non-threatened",
    red_list_category_code %in% c("EW")                                        ~ "EW",
    red_list_category_code %in% c("Ex?", "EX", "E", "Ex", "Ex/E")              ~ "EX",
    red_list_category_code %in% c("DD")                                        ~ "DD",
    red_list_category_code %in% c("NE", "N/A")                                 ~ NA, #"Not evaluated",
    red_list_category_code %in% c("T", "R", "I", "K")                          ~ NA  #"Unknown"
  )) %>%
  mutate(year_published = case_when(is.na(red_list_category_code) ~ NA,
                                    .default = year_published)) %>%
  filter(!(is.na(year_published) & is.na(red_list_category_code)))

### List of non-extant species - some species are extinct in old assessments but not in the latest assessment
extinctSp <- all_assessments %>%
  group_by(taxon_scientific_name) %>%
  filter(any(red_list_category_code %in% c("EW", "EX"))) %>%
  filter(all(red_list_category_code %in% c(NA, "EW", "EX")) |
           red_list_category_code %in% c("EW", "EX") & latest == TRUE) %>%
  pull(taxon_scientific_name) %>%
  unique()

### Remove extinct species - 65,071 assessments for 36,857 species
all_assessments <- all_assessments %>%
  filter(!taxon_scientific_name %in% extinctSp)
table(all_assessments$red_list_category_code, useNA = "always")
nrow(all_assessments)
length(unique(all_assessments$taxon_scientific_name))

### For ordering columns later
assessment_years <- c("latest_RL_cat", "latest_RL_year",
                      sort(unique(all_assessments$year_published)))

### Get most recent RL status - 36,857 species
latest_assessments <- all_assessments %>%
  arrange(taxon_scientific_name, -year_published) %>%
  group_by(taxon_scientific_name) %>%
  mutate(latest_RL_cat  = case_when(all( is.na(year_published)) ~ "Not evaluated",
                                    any(!is.na(year_published)) ~ red_list_category_code[1])) %>%
  mutate(latest_RL_year = case_when(all( is.na(year_published)) ~ NA,
                                    any(!is.na(year_published)) ~ year_published[1])) %>%
  select(taxon_scientific_name, latest_RL_cat, latest_RL_year) %>%
  distinct()
table(latest_assessments$latest_RL_cat, useNA = "always")

### Pivot assessment info wider
assessments_wide <- all_assessments %>%
  arrange(taxon_scientific_name, -year_published) %>%
  group_by(taxon_scientific_name) %>%
  pivot_wider(id_cols     = taxon_scientific_name,
              names_from  = year_published,
              values_from = red_list_category_code,
              values_fill = NA)

### Join together with latest assessments
assessments_wide <- left_join(latest_assessments, assessments_wide, by = "taxon_scientific_name") %>%
  select(taxon_scientific_name, all_of(assessment_years)) %>%
  mutate(latest_RL_cat = case_when(is.na(latest_RL_cat) ~ "Not evaluated",
                                   .default = as.factor(latest_RL_cat)))

### Join to distribution info
threatLevels <- left_join(all_intersections,
                          assessments_wide,
                          join_by("NameIUCN" == "taxon_scientific_name")) %>%
  rename(Species = NameIUCN)

### Which species are subregion endemics
threatLevels <- threatLevels %>%
  mutate(subregionEndemic = case_when(AsiaEndemic == 1 & rowSums(select(., Indian_Subcontinent:New_Guinea)) == 1 ~ 1,
                                      AsiaEndemic == 1 & rowSums(select(., Andaman:South_Kuroshio))         == 1 ~ 1,
                                      .default = 0)) %>%
  mutate(latest_RL_cat = case_when(is.na(latest_RL_cat) ~ "Not evaluated",
                                   .default = latest_RL_cat)) %>%
  mutate(latest_RL_cat = factor(latest_RL_cat,
                                levels = c("CR", "EN", "VU", "Non-threatened", "DD", "Not evaluated"))) %>%
  mutate(Group = factor(Group,
                        levels = rev(unique(taxon_info$Group)))) %>%
  mutate(TaxonForPlot = factor(TaxonForPlot,
                               levels = rev(unique(taxon_info$TaxonForPlot))))

table(threatLevels$latest_RL_cat,  useNA = "always")
#    CR             EN             VU Non-threatened             DD  Not evaluated           <NA> 
#  1943           3780           3463          22444           5698         110396              0

table(threatLevels$latest_RL_year, useNA = "always")
# 1996   1998   2000   2004   2006   2007   2008   2009   2010   2011   2012   2013   2014   2015
#   68    724     38    108      1      4    683     64    907    832   1240    520    626    651 
# 2016   2017   2018   2019   2020   2021   2022   2023   2024   2025   2026   <NA>
# 1546    811   2397   4000   5077   4106   2769   2116   6177   1856      7 110396

### Prop. species assessed / without assessments
reframe(threatLevels,
        n = n(),
        n_not_assessed = sum(latest_RL_cat == "Not evaluated"),
        n_assessed     = sum(latest_RL_cat != "Not evaluated"),
        perc_asessed   = (n_assessed / n) * 100,
        .by = TaxonForPlot) %>%
  print(n = 23)
#    TaxonForPlot               n n_not_assessed n_assessed perc_asessed
# 1  Ferns                   3903           3736        167      4.28   
# 2  Flowering plants       70082          53011      17071     24.4    
# 3  Gymnosperms              288             47        241     83.7    
# 4  Lycophytes               389            387          2      0.514  
# 5  Amphibians              1885              0       1885    100      
# 6  Birds                   2918              0       2918    100      
# 7  Freshwater fish         3205              0       3205    100      
# 8  Mammals                 1501              1       1500     99.9    
# 9  Reptiles                2500              1       2499    100.0    
# 10 Ants                    2971           2968          3      0.101  
# 11 Bees                    2134           2131          3      0.141  
# 12 Butterflies             4757           4529        228      4.79   
# 13 Caddisflies             5873           5873          0      0      
# 14 Centipedes               488            488          0      0      
# 15 Flies                  24414          24413          1      0.00410
# 16 Freshwater crabs         670              0        670    100      
# 17 Millipedes              2287           2282          5      0.219  
# 18 Mirid bugs              1440           1440          0      0      
# 19 Spiders                 7624           7594         30      0.393  
# 20 Stick and leaf insects  1495           1494          1      0.0669 
# 21 Bony fish               5789              1       5788    100.0    
# 22 Reef corals              708              0        708    100      
# 23 Sharks and rays          403              0        403    100  

### By group
reframe(threatLevels,
        n = n(),
        n_not_assessed = sum(latest_RL_cat == "Not evaluated"),
        n_assessed     = sum(latest_RL_cat != "Not evaluated"),
        perc_asessed   = (n_assessed / n) * 100,
        .by = GroupForPlot)
#   GroupForPlot                                   n n_not_assessed n_assessed perc_asessed
# 1 "Vascular\nplants"                         74662          57181      17481        23.4 
# 2 "Terrestrial\n& freshwater\nvertebrates"   12009              2      12007       100.0 
# 3 "Terrestrial\n& freshwater\ninvertebrates" 54153          53212        941         1.74
# 4 "Marine"                                    6900              1       6899       100.0

### Proportion assessed and threatened for region as a whole
reframe(threatLevels, No_Threat = n(), .by = c(latest_RL_cat)) %>%
  mutate(Prop_threat = No_Threat / sum(No_Threat))
#   latest_RL_cat  No_Threat Prop_threat
# 1 Not evaluated     110396      0.747 
# 2 Non-threatened     22444      0.152 
# 3 DD                  5698      0.0386
# 4 CR                  1943      0.0132
# 5 EN                  3780      0.0256
# 6 VU                  3463      0.0234

### Proportion assessed and threatened for region as a whole
iucn_summary <- threatLevels %>%
  reframe(No_Threat = n(),
          .by = c(latest_RL_cat, TaxonForPlot)) %>%
  pivot_wider(id_cols     = TaxonForPlot,
              names_from  = latest_RL_cat,
              values_from = No_Threat) %>%
  rename(Taxon = TaxonForPlot) %>%
  select(Taxon, CR, EN, VU, `Non-threatened`, DD, `Not evaluated`)
iucn_summary[is.na(iucn_summary)] <- 0

write_csv(iucn_summary, file.path(threatDir, "iucn_summary.csv"))

####################################################################################################
### FIGURE 5A: Prop threatened for tropical Asia
####################################################################################################

group_colours <- taxon_info %>%
  reframe(col = rgb(apply(col2rgb(Colour), 1, mean)[1] / 255,
                    apply(col2rgb(Colour), 1, mean)[2] / 255,
                    apply(col2rgb(Colour), 1, mean)[3] / 255),
          .by = Group)

### Perc. species that are endemics
percSp <- threatLevels %>%
  reframe(perc = paste0(format(round(sum(subregionEndemic == 1) / n() * 100, 1), nsmall = 1), "%"),
          .by = c(TaxonForPlot))

theme_propRL <- theme(panel.background  = element_blank(),
                      panel.border      = element_blank(),
                      panel.grid        = element_blank(),
                      plot.margin       = margin(t = 0, r = 0, b = 0, l = 5),
                      axis.ticks        = element_line(colour = "black", linewidth = 0.25),
                      axis.ticks.length = unit(0.5, units = "mm"),
                      axis.text         = element_text(family = "Arial", colour = "black", size = 5),
                      axis.title        = element_text(family = "Arial", colour = "black", size = 5),
                      title             = element_text(family = "Arial", colour = "black", size = 5))

### All
p1 <- threatLevels %>%
  filter(subregionEndemic == 1) %>%
  ggplot(aes(x = TaxonForPlot, fill = latest_RL_cat)) +
  theme_propRL +
  coord_flip(clip = "off") +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.75)),
                   limits = c(rev(taxon_info$TaxonForPlot[taxon_info$Group == "Marine"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Terrestrial and freshwater invertebrates"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Terrestrial and freshwater vertebrates"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Vascular plants"]))) +
  scale_y_continuous(expand = expansion(add = c(0, 0))) +
  labs(y = "Prop. sp.", x = NULL, title = "Species found in\nonly one subregion") +
  annotate("rect", xmin = 22, xmax = Inf, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Vascular plants"]) +
  annotate("rect", xmin = 16, xmax = 22, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Terrestrial and freshwater vertebrates"]) +
  annotate("rect", xmin = 4, xmax = 16, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Terrestrial and freshwater invertebrates"]) +
  annotate("rect", xmin = -Inf, xmax = 4, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Marine"]) +
  scale_fill_manual(values = colsRedList) +
  geom_bar(stat = "count", position = position_fill(reverse = FALSE)) +
  annotate("text", label = percSp$perc, x = percSp$TaxonForPlot, y = 1,
           colour = "black", family = "Arial", size.unit = "pt", size = 5, hjust = 1) +
  annotate("rect", xmin = -0.01, xmax = 27.01, ymin = 0, ymax = 1,
           fill = NA, colour = "black", linewidth = 0.25)

p2 <- threatLevels %>%
  filter(subregionEndemic == 0) %>%
  ggplot(aes(x = TaxonForPlot, fill = latest_RL_cat)) +
  theme_propRL +
  coord_flip(clip = "off") +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.75)),
                   limits = c(rev(taxon_info$TaxonForPlot[taxon_info$Group == "Marine"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Terrestrial and freshwater invertebrates"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Terrestrial and freshwater vertebrates"]), "",
                              rev(taxon_info$TaxonForPlot[taxon_info$Group == "Vascular plants"]))) +
  scale_y_continuous(expand = expansion(add = c(0, 0))) +
  labs(y = "Prop. sp.", x = NULL, title = "All other species") +
  annotate("rect", xmin = 22, xmax = Inf, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Vascular plants"]) +
  annotate("rect", xmin = 16, xmax = 22, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Terrestrial and freshwater vertebrates"]) +
  annotate("rect", xmin = 4, xmax = 16, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Terrestrial and freshwater invertebrates"]) +
  annotate("rect", xmin = -Inf, xmax = 4, ymin = 0, ymax = 1, alpha = 0.5, colour = NA,
           fill = group_colours$col[group_colours$Group == "Marine"]) +
  scale_fill_manual(values = colsRedList) +
  geom_bar(stat = "count", position = position_fill(reverse = FALSE)) +
  annotate("rect", xmin = -0.01, xmax = 27.01, ymin = 0, ymax = 1,
           fill = NA, colour = "black", linewidth = 0.25)

### Subregion endemics and non-endemics only
propRL <- p1 + p2 +
  plot_layout(ncol = 2, guides = "collect", axes = "collect") &
  theme(legend.position = "none")

ggsave(propRL,
       filename = file.path(figDir, "prop_threatened_endemics.png"),
       width = 60, height = 80, units = "mm", dpi = 600, bg = "white")

ggsave(propRL,
       filename = file.path(figDir, "prop_threatened_endemics.pdf"),
       width = 60, height = 80, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### FIGURE 5B: Prop threatened by subregion - map version - All taxa
####################################################################################################

### Get subregion centroids
coords <- bind_rows(st_drop_geometry(terr_bioregions),
                    rename(st_drop_geometry(mare_bioregions), Bioregion = PROVINCE)) %>%
  mutate(Bioregion = gsub(" ", "_", Bioregion)) %>%
  select(Bioregion, X, Y) %>%
  mutate(X = case_when(Bioregion == "Sulawesi"                     ~ 13490000,
                       Bioregion == "Andaman"                      ~ 10400000,
                       Bioregion == "Bay_of_Bengal"                ~  9600000,
                       Bioregion == "Central_Indian_Ocean_Islands" ~  8150000,
                       Bioregion == "Eastern_Coral_Triangle"       ~ 16700000,
                       Bioregion == "Sahul_Shelf"                  ~ 15100000,
                       Bioregion == "South_China_Sea"              ~ 12800000,
                       Bioregion == "West_and_South_Indian_Shelf"  ~  9000000,
                       .default = X)) %>%
  mutate(Y = case_when(Bioregion == "Sumatra"                     ~ -370000,
                       Bioregion == "Sulawesi"                    ~ -180000,
                       Bioregion == "Andaman"                     ~ 1100000,
                       Bioregion == "Eastern_Coral_Triangle"      ~ -600000,
                       Bioregion == "South_China_Sea"             ~ 1750000,
                       Bioregion == "West_and_South_Indian_Shelf" ~  760000,
                       .default = Y))

### Get tropical Asia richness for standardisation
regRich <- threatLevels %>%
  reframe(RichRegion     = n(),
          RichRegionEval = sum(latest_RL_cat != "Not evaluated"),
          .by = GroupForPlot)

### Get subregions richness for standardisation
subregRich <- threatLevels %>%
  reframe(across(Indian_Subcontinent:South_Kuroshio, \(x) sum(x, na.rm = TRUE)), .by = GroupForPlot)  %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio,
               names_to = "Subregion",
               values_to = "RichSubregion")

########################################################
### Prop. threatened - Tropical Asia as a whole

### Summarise each taxon for tropical Asia - all species
propThreatAll <- threatLevels %>%
  select(GroupForPlot, Species, latest_RL_cat) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, latest_RL_cat)) %>%
  mutate(redlistCategory = factor(latest_RL_cat,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(regRich, by = "GroupForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(LogRichRegion = log(RichRegion)) %>%
  mutate(Subset = "All species")

### Summarise each taxon for tropical Asia - only species with IUCN evaluations
propThreatEval <- threatLevels %>%
  filter(latest_RL_cat != "Not evaluated") %>%
  select(GroupForPlot, Species, latest_RL_cat) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, latest_RL_cat)) %>%
  mutate(redlistCategory = factor(latest_RL_cat,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(regRich, by = "GroupForPlot") %>%
  mutate(propThreat    = NoThreat / RichRegionEval) %>%
  mutate(LogRichRegion = log(RichRegionEval)) %>%
  mutate(Subset = "Evaluated\nspecies only")

### Combine
propThreat <- bind_rows(propThreatAll, propThreatEval) %>%
  mutate(Group = paste(GroupForPlot, Subset))

### Loop through and plot pie charts for tropical Asia
propThreat_group_plot_list <- list()
for(i in 1:length(unique(propThreat$Group))) {
  ### Subset taxon
  taxon <- unique(propThreat$Group)[i]
  dd <- filter(propThreat, Group == taxon)
  
  ### Generate plot
  p <- ggplot(dd,
              aes(x = LogRichRegion / 2,
                  y = NoThreat,
                  fill = redlistCategory,
                  width = 11.1
              )) +
    theme_void() +
    facet_grid(Subset ~ GroupForPlot, switch = "y") +
    geom_col() +
    scale_fill_manual(values = colsRedList) +
    coord_polar("y", start = 0)
  
  ### Taxonomic group
  if(i %in% c(1:4)) {
    p <- p +
      theme(strip.text.x = element_text(family = "arial", colour = "black", size = 8))
  } else {
    p <- p +
      theme(strip.text.x = element_blank())
  }
  
  ### Species subset
  if(i %in% c(1, 5)) {
    p <- p +
      theme(strip.text.y = element_text(family = "arial", colour = "black", size = 8, angle = 90))
  } else {
    p <- p +
      theme(strip.text.y = element_blank())
  }
  if(i %in% c(1:4)) { p <- p + labs(x = "All",   y = "All",   title = "All") }
  if(i %in% c(5:8)) { p <- p + labs(x = "Eval.", y = "Eval.", title = "Eval.") }
  
  ### Legend
  if(i == 1) {
    p <- p +
      theme(legend.text        = element_text(family = "arial", colour = "black", size = 6),
            legend.title       = element_text(family = "arial", colour = "black", size = 8),
            legend.key.spacing = unit(3, units = "mm")) +
      guides(fill = guide_legend(title = "Red List category", position = "bottom",
                                 nrow = 1, reverse = TRUE))
  } else {
    p <- p +
      guides(fill = guide_none())
  }
  propThreat_group_plot_list[[i]] <- p
}

### Name list for plotting on top of maps later on
names(propThreat_group_plot_list) <- c("All_plants", "All_verts", "All_inverts", "All_marine",
                                       "Eval_plants", "Eval_verts", "Eval_inverts", "Eval_marine")

########################################################
### Prop. threatened - by subregion

### Summarise across subregions
subregionIucn <- threatLevels %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio, names_to = "Subregion", values_to = "Presence") %>%
  relocate(Subregion) %>%
  filter(Presence == 1) %>%
  select(GroupForPlot, Subregion, latest_RL_cat) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, Subregion, latest_RL_cat)) %>%
  mutate(redlistCategory = factor(latest_RL_cat,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(subregRich, by = c("GroupForPlot", "Subregion")) %>%
  mutate(propThreat = NoThreat / RichSubregion) %>%
  mutate(LogRichSubregion = log(RichSubregion)) %>%
  left_join(select(region_info, Bioregion), join_by("Subregion" == "Bioregion")) %>%
  left_join(coords, join_by("Subregion" == "Bioregion")) %>%
  pivot_wider(id_cols = c(GroupForPlot, Subregion, RichSubregion, LogRichSubregion, X, Y),
              names_from = redlistCategory,
              values_from = propThreat)
subregionIucn[is.na(subregionIucn)] <- 0

### Proportion non-evaluated species in subregions
subregionIucn %>%
  mutate(n_not_assessed = RichSubregion * `Not evaluated`) %>%
  mutate(n_dd           = RichSubregion * DD) %>%
  reframe(RichSubregion        = sum(RichSubregion, na.rm = TRUE),
          n_not_assessed       = sum(n_not_assessed, na.rm = TRUE),
          n_dd                 = sum(n_dd, na.rm = TRUE),
          perc_not_assessed    = (n_not_assessed / RichSubregion) * 100,
          perc_dd              = (n_dd / RichSubregion) * 100,
          perc_dd_not_assessed = ((n_dd + n_not_assessed) / RichSubregion) * 100,
          .by = Subregion) %>%
  print(n = 23)

### Proportion threatened species in subregions
threatLevels %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio, names_to = "Subregion", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  select(GroupForPlot, Subregion, latest_RL_cat) %>%
  reframe(no_assessed = n(),
          .by = c(GroupForPlot, Subregion, latest_RL_cat)) %>%
  left_join(subregRich, by = c("GroupForPlot", "Subregion")) %>%
  left_join(select(region_info, Bioregion), join_by("Subregion" == "Bioregion")) %>%
  pivot_wider(id_cols = c(GroupForPlot, Subregion, RichSubregion),
              names_from = latest_RL_cat,
              values_from = no_assessed) %>%
  mutate(non_threatened = `Non-threatened` + DD) %>%
  mutate(threatened     = CR + EN + VU) %>%
  reframe(non_threatened = sum(non_threatened, na.rm = TRUE),
          threatened     = sum(threatened, na.rm = TRUE),
          .by = Subregion) %>%
  mutate(perc_threatened = (threatened / non_threatened) * 100) %>%
  print(n = 23)

### Loop through taxonomic groups and make maps
facet_groups <- unique(subregionIucn$GroupForPlot)
region_iucn_plot_list <- list()
for(i in 1:4) {
  ### Get data for group
  dd <- filter(subregionIucn, GroupForPlot == facet_groups[i])
  
  ### Terrestrial groups
  if(i %in% 1:3) {
    basePlot <- ggplot() + 
      theme(panel.grid = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm"),
            panel.background = element_rect(fill = "white", colour = NA)) +
      geom_sf(data = base_map, color = NA, fill = "grey50") +
      geom_sf(data = terr_bioregions, color = NA, fill = "grey15") +
      geom_sf(data = terr_divides, color = "grey90", linewidth = 0.2) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                col = "black", fill = NA) +
      geom_text(aes(x = 7800000, y = -1600000),#x = 7714253, y = -1688714),
                label = facet_groups[i],
                hjust = 0, vjust = 0,
                color = "black", size.unit = "pt", size = 7, family = "sans", fontface = "plain")
  }
  
  ### Marine groups
  if(i == 4) {
    basePlot <- ggplot() + 
      theme(panel.grid = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm"),
            panel.background = element_rect(fill = "grey50", colour = NA)) +
      geom_sf(data = base_map, color = NA, fill = "white") +
      geom_sf(data = mare_bioregions, color = NA, fill = "grey15") +
      geom_sf(data = mare_divides, color = "grey90", linewidth = 0.2) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                col = "black", fill = NA) +
      geom_text(aes(x = 7800000, y = -1600000),#x = 7714253, y = -1688714),
                label = facet_groups[i],
                hjust = 0, vjust = 0,
                color = "white", size.unit = "pt", size = 7, family = "sans", fontface = "plain")
  }
  
  ### Add in base pie chart
  region_iucn_plot_list[[i]] <- basePlot
  
  ### Add in pie charts
  region_iucn_plot_list[[i]] <-  region_iucn_plot_list[[i]] +
    geom_scatterpie(data = dd,
                    aes(x = X, y = Y, group = Subregion,
                        # r = LogRichSubregion * 38000
                        r = max(subregionIucn$LogRichSubregion) * 38000
                    ),
                    cols = c("CR", "EN", "VU", "Non-threatened", "DD", "Not evaluated"),
                    colour = NA) +
    scale_fill_manual(values = colsRedList,
                      name = "Red List category") +
    coord_sf(xlim = c(7614253, 17232260), ylim = c(-1788714, 3279685))
  
  ### And RL legend
  if(i == 1) {
    region_iucn_plot_list[[i]] <- region_iucn_plot_list[[i]] +
      guides(fill = guide_legend(title = "Red List category", position = "bottom",
                                 nrow = 1, reverse = FALSE))
  } else {
    region_iucn_plot_list[[i]] <- region_iucn_plot_list[[i]] +
      guides(fill = guide_none())
  }
}

### Add in subplots (prop. threatened across all tropical Asia)
pBox <- ggplot() + 
  theme(plot.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.background = element_blank())
pBoxMarine <- ggplot() + 
  theme(plot.background = element_rect(fill = "grey50", colour = "black", linewidth = 0.5),
        panel.background = element_blank())

pSubPlants <-
  (propThreat_group_plot_list[["All_plants"]]  + guides(fill = "none")) +
  (propThreat_group_plot_list[["Eval_plants"]] + guides(fill = "none")) &
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5),
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[1]] <- region_iucn_plot_list[[1]] + 
  inset_element(pBox,       left = 0.71, right = 0.98, bottom = 0.60, top = 0.94) +
  inset_element(pSubPlants, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubVerts <- 
  propThreat_group_plot_list[["All_verts"]] + propThreat_group_plot_list[["Eval_verts"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[2]] <- region_iucn_plot_list[[2]] + 
  inset_element(pBox,      left = 0.71, right = 0.98, bottom = 0.60, top = 0.94) +
  inset_element(pSubVerts, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubInverts <- 
  propThreat_group_plot_list[["All_inverts"]] + propThreat_group_plot_list[["Eval_inverts"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[3]] <- region_iucn_plot_list[[3]] + 
  inset_element(pBox,        left = 0.71, right = 0.98, bottom = 0.60, top = 0.94) +
  inset_element(pSubInverts, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubMarine <-
  propThreat_group_plot_list[["All_marine"]] + propThreat_group_plot_list[["Eval_marine"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[4]] <- region_iucn_plot_list[[4]] +
  inset_element(pBoxMarine, left = 0.72, right = 0.99, bottom = 0.46, top = 0.79) +
  inset_element(pSubMarine, left = 0.72, right = 0.99, bottom = 0.45,  top = 0.84)

pIucnCatsSubregionMap <-
  (region_iucn_plot_list[[1]] +  region_iucn_plot_list[[2]]) /
  (region_iucn_plot_list[[3]] +  region_iucn_plot_list[[4]]) +
  plot_layout(tag_level = "new", guides = "collect") &
  theme(legend.position = "bottom",
        legend.title    = element_text(family = "Arial", colour = "black", size = 6),
        legend.text     = element_text(family = "Arial", colour = "black", size = 5))

ggsave(pIucnCatsSubregionMap,
       filename = file.path(figDir, "map_prop_threatened.png"),
       width = 120, height = 80, units = "mm", dpi = 600, bg = "white")

ggsave(pIucnCatsSubregionMap,
       filename = file.path(figDir, "map_prop_threatened.pdf"),
       width = 120, height = 80, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### FIGURE 5C: HFP against latest assessment date
####################################################################################################

### Group - subregion richness
subregionRichGroup <- threatLevels %>%
  pivot_longer(cols      = Indian_Subcontinent:New_Guinea,
               names_to  = "Bioregion",
               values_to = "presence") %>%
  filter(presence == 1) %>%
  reframe(nsp = length(unique(Species)),
          .by = c(Group, GroupForPlot, Bioregion))

### Taxon - subregion no. assessments per year
subregionAssessGroup <- threatLevels %>%
  pivot_longer(cols      = Indian_Subcontinent:New_Guinea,
               names_to  = "Bioregion",
               values_to = "presence") %>%
  filter(presence == 1) %>%
  reframe(n_assess_group = sum(!is.na(latest_RL_year)),
          .by = c(GroupForPlot, Bioregion, latest_RL_year)) %>%
  filter(!is.na(latest_RL_year)) %>%
  rename(Year = latest_RL_year) %>%
  arrange(GroupForPlot, Bioregion, Year)

### Add in missing years with no taxon-bioregion assessments
subregionAssessGroup <- full_join(subregionAssessGroup,
                                  expand_grid(GroupForPlot = unique(subregionAssessGroup$GroupForPlot),
                                              Bioregion    = unique(subregionAssessGroup$Bioregion),
                                              Year         = min(subregionAssessGroup$Year):max(subregionAssessGroup$Year)),
                                  by = c("GroupForPlot", "Bioregion", "Year")) %>%
  arrange(GroupForPlot, Bioregion, Year) %>%
  mutate(n_assess = case_when(is.na(n_assess_group) ~ 0, .default = n_assess_group))

### Add in subregion richness and calculate cumulative proportion species described
subregionAssessGroup <- full_join(subregionAssessGroup,
                                  subregionRichGroup,
                                  by = c("GroupForPlot", "Bioregion")) %>%
  group_by(GroupForPlot, Bioregion) %>%
  mutate(prop_assess_group = n_assess_group / nsp) %>%
  ungroup()

### Merge with HFP and standardise for use on secondary axis
maxProps <- reframe(subregionAssessGroup,
                    max_prop = max(prop_assess_group, na.rm = TRUE),
                    .by = GroupForPlot)
maxProps <- c("Terrestrial\n& freshwater\ninvertebrates" = 0.0181,
              "Terrestrial\n& freshwater\nvertebrates"   = 0.314, 
              "Vascular\nplants"                         = 0.105)

subregionAssessGroup <- subregionAssessGroup %>%
  left_join(hfp, by = c("Bioregion", "Year")) %>%
  left_join(maxProps, by = "GroupForPlot") %>%
  mutate(mean_hfp_stand = (mean_hfp * max_prop) / 35,
         q25_hfp_stand  = (q25_hfp  * max_prop) / 35,
         q75_hfp_stand  = (q75_hfp  * max_prop) / 35)

### Add in bioregion data and tidy up factors
subregionAssessGroup <- subregionAssessGroup %>%
  select(-Realm, -LabelForPlot, -Label, -Order, -Colour) %>%
  left_join(region_info, by = "Bioregion") %>%
  mutate(LabelForPlot = factor(LabelForPlot,
                               levels = region_info$LabelForPlot)) %>%
  mutate(GroupForPlot = factor(GroupForPlot,
                               levels = unique(taxon_info$GroupForPlot)))

### Plot theme
theme_assess_year <- theme(
  plot.margin        = margin(t = 0, r = 0, b = 2, l = 0, "mm"),
  panel.spacing      = unit(1, "mm"),
  panel.background   = element_blank(),
  panel.border       = element_rect(fill = NULL, colour = "black", linewidth = 0.5),
  panel.grid         = element_blank(),
  axis.ticks         = element_line(colour = "black", linewidth = 0.25),
  axis.ticks.length  = unit(0.5, units = "mm"),
  axis.text.x        = element_text(family = "Arial", colour = "black", size = 5, angle = 90, vjust = 0.5),
  axis.text.y        = element_text(family = "Arial", colour = "black", size = 5),
  axis.title.x       = element_text(family = "Arial", colour = "black", size = 6),
  axis.title.y       = element_text(family = "Arial", colour = "black", size = 6,
                                    margin = unit(c(t = 0, r = 0, b = 0, l = 0), "mm")),
  strip.background.x = element_rect(fill = "grey30", colour = "black", linewidth = 0.25),
  strip.text.x       = element_text(family = "Arial", colour = "white", size = 6)
) 

p_list <- lapply(
  split(subregionAssessGroup, subregionAssessGroup$GroupForPlot),
  FUN = function(x) {
    p <- ggplot(data = x) +
      theme_assess_year +
      facet_grid( ~ LabelForPlot, switch = "y") +
      labs(x = "Year", y = "Prop. species with latest assessements") +
      scale_x_continuous(breaks = seq(2000, 2024, 8)) +
      scale_y_continuous(sec.axis = sec_axis(name = "Human Footprint",
                                             transform = ~ . * (35 / max(x$prop_assess_group, na.rm = TRUE))
      )) +
      geom_polygon(data = reframe(filter(x, !is.na(q25_hfp_stand)),
                                  year = c(Year, rev(Year)),
                                  iqr  = c(q25_hfp_stand, rev(q75_hfp_stand)),
                                  .by  = c(GroupForPlot, LabelForPlot)),
                   aes(year, iqr),
                   fill = "grey80") +
      geom_line(aes(Year, mean_hfp_stand), linewidth = 0.3) +
      geom_bar(aes(x = Year, y = prop_assess_group), stat = "identity",
               width = 1, fill = "red4", alpha = 1)
    return(p)
  }
)
plotLatestAssess <- p_list[[1]] +
  (p_list[[2]] + theme(strip.background = element_blank(), strip.text = element_blank())) +
  (p_list[[3]] + theme(strip.background = element_blank(), strip.text = element_blank())) +
  plot_layout(nrow = 3, axes = "collect", axis_titles = "collect")

####################################################################################################
### EXTENDED DATA FIGURE 4: Prop threatened for tropical Asia - each taxon separately
####################################################################################################

### Get regional richness for pie size
regRich <- threatLevels %>%
  reframe(RichRegion = n(), .by = TaxonForPlot)

### Summarise each taxon
propThreat <- threatLevels %>%
  select(TaxonForPlot, GroupForPlot, Species, latest_RL_cat) %>%
  reframe(NoThreat = n(),
          .by = c(TaxonForPlot, GroupForPlot, latest_RL_cat)) %>%
  mutate(latest_RL_cat = factor(latest_RL_cat,
                                levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR", "EX"))) %>%
  left_join(regRich, by = "TaxonForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(LogRichRegion = log(RichRegion))

for(i in 1:length(unique(propThreat$TaxonForPlot))) {
  ### Subset taxon
  taxon <- unique(threatLevels$TaxonForPlot)[i]
  dd <- filter(propThreat, TaxonForPlot == taxon)
  
  ### Generate plot
  p <- ggplot(dd,
              aes(x = LogRichRegion / 2,
                  y = NoThreat,
                  fill = latest_RL_cat,
                  width = 11.1
              )) +
    theme_void() +
    theme(strip.text.y   = element_text(family = "arial", colour = "black", size = 8, angle = 90),
          strip.text.x   = element_text(family = "arial", colour = "black", size = 8),
          legend.text    = element_text(family = "arial", colour = "black", size = 6),
          legend.title   = element_text(family = "arial", colour = "black", size = 8),
          legend.key.spacing = unit(3, units = "mm")
    ) +
    facet_wrap( ~ TaxonForPlot) +
    geom_col() +
    scale_fill_manual(values = colsRedList) +
    coord_polar("y", start = 0)
  if(i %in% c(1, 5, 10, 16, 21)) {
    p <- p +
      facet_grid(GroupForPlot ~ TaxonForPlot, switch = "y")
  }
  if(i == 1) {
    p <- p +
      guides(fill = guide_legend(title = "Red List category", position = "bottom",
                                 nrow = 1, reverse = TRUE))
  } else {
    p <- p +
      guides(fill = guide_none())
  }
  assign(paste0("p", i), p)
}

pIucnCatsAll <- p1 + p2 + p3 + p4 + plot_spacer() + plot_spacer() +
  p5 + p6 + p7 + p8 + p9 + plot_spacer() +
  p10 + p11 + p12 + p13 + p14 + p15 +
  p16 + p17 +  p18 + p19 + p20 + plot_spacer() +
  p21 + p22 + p23 + 
  plot_layout(guides = "collect", tag_level = "keep", nrow = 5, ncol = 6, axis_titles = "collect") &
  theme(legend.position = "bottom")
# pIucnCatsAll

ggsave(file.path(figDir, "Prop_threatened.png"), pIucnCatsAll,
       width = 180, height = 175, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "Prop_threatened.pdf"), pIucnCatsAll,
       width = 180, height = 175, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### EXTENDED DATA FIGURE 5: Year published - each taxon separately
####################################################################################################

for(i in 1:length(unique(threatLevels$TaxonForPlot))) {
  ### Subset taxon
  taxon <- unique(threatLevels$TaxonForPlot)[i]
  dd <- droplevels(filter(threatLevels, TaxonForPlot == taxon)) %>%
    filter(!is.na(latest_RL_year))
  
  ### If there are some IUCN data for taxon
  if(nrow(dd) > 0) {
    p <- ggplot(dd, aes(latest_RL_year)) +
      theme(panel.background = element_rect(fill = NA, colour = "black"),
            axis.text   = element_text(family = "Arial", size = 6, colour = "black"),
            axis.title  = element_text(family = "Arial", size = 8, colour = "black"),
            axis.ticks  = element_line(colour = "black", linewidth = 0.25),
            axis.line   = element_line(colour = "black", linewidth = 0.25),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey50", linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = "grey50", linewidth = 0.1)) +
      geom_text(data = data.frame(x = 1996, y = max(table(dd$latest_RL_year)) * 0.95,
                                  label = paste0(unique(dd$TaxonForPlot), "\n",
                                                 round((sum(dd$latest_RL_year < 2015) / nrow(dd)) * 100, 1),
                                                 "%")),
                aes(x = x, y = y, label = label),vjust = "inward", hjust = "inward",
                family = "Arial", colour = "black", size = 6, size.unit = "pt") +
      scale_x_continuous(breaks = function(x) unique(floor(pretty(x, 3))),
                         expand = expansion(mult = c(0.02, 0.02)), limits = c(1996, 2025)) +
      scale_y_continuous(breaks = function(x) unique(floor(pretty(x, 4))),
                         expand = expansion(mult = c(0, 0.05)),
                         labels = label_number(accuracy = 1, big.mark = "")) +
      geom_vline(xintercept = 2015, linetype = 2, col = "red") +
      geom_histogram(breaks = seq(1995.5, 2025.5, 1))
  }
  
  ### If there is no IUCN data generate blank plot
  if(nrow(dd) == 0) {
    dd <- data.frame(TaxonForPlot = taxon)
    p <- ggplot(dd) +
      theme(panel.background = element_rect(fill = NA, colour = "black"),
            axis.text   = element_text(family = "Arial", size = 6, colour = "black"),
            axis.title  = element_text(family = "Arial", size = 8, colour = "black"),
            axis.ticks  = element_line(colour = "black", linewidth = 0.25),
            axis.line   = element_line(colour = "black", linewidth = 0.25),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey50", linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = "grey50", linewidth = 0.1)) +
      geom_text(data = data.frame(x = 1996, y = 1 * 0.95,
                                  label = unique(dd$TaxonForPlot)),
                aes(x = x, y = y, label = label),vjust = "inward", hjust = "inward",
                family = "Arial", colour = "black", size = 6, size.unit = "pt") +
      scale_y_continuous(labels = NULL, breaks = NULL, limits = c(0, 1)) +
      lims(x = c(1995.5, 2025.5))
  }
  
  ### Add axis labels to specific plots
  if(i %in% c(1, 5, 10, 16, 21)) {
    p <- p
  }
  if(i %in% c(15)) {
    p <- p + labs(x = "Assessment year")
  } else {
    p <- p + labs(x = NULL)
  }
  if(i %in% c(23)) {
    p <- p + labs(y = "No. of species")
  } else {
    p <- p + labs(y = NULL)
  }
  assign(paste0("p", i), p)
}
pAssessmentYear <-
  p1 +            p5 +            p10 +            p16 +           p21 +
  p2 +            p6 +            p11 +            p17 +           p22 +
  p3 +            p7 +            p12 +            p18 +           p23 +
  p4 +            p8 +            p13 +            p19 +           plot_spacer() +
  plot_spacer() + p9 +            p14 +            p20 +           plot_spacer() +
  plot_spacer() + plot_spacer() + p15 +            plot_spacer() + plot_spacer() +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect_x", tag_level = "keep",
              nrow = 6, ncol = 5)
# pAssessmentYear

ggsave(file.path(figDir, "Assessment_year_by_taxon.png"), pAssessmentYear,
       width = 180, height = 135, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "Assessment_year_by_taxon.pdf"), pAssessmentYear,
       width = 220, height = 135, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)
