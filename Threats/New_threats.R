####################################################################################################
### 
### Charlie Marsh
### charliem2003@github
### 05/2025
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
library(scales)
library(ggplot2)
library(scatterpie)
library(patchwork)
library(sf)

### Locations of data, scripts and results
projDir  <- file.path("/mnt", "Work", "NUS", "BTAS")
interDir <- file.path(projDir, "Threats")
mapDir   <- file.path(projDir, "Maps")
figDir   <- file.path(projDir, "Threats", "Figures")


### Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
projDir  <- "Threats"                                             # project dir

### You shouldn't need to adjust these folders
mapDir   <- file.path(projDir, "Data")
figDir   <- file.path(projDir, "Figures")
if(!dir.exists(figDir)) { dir.create(figDir) }

#==================================================================================================#
#----------------------------------------- Read in data -------------------------------------------#
#==================================================================================================#

### Combined intersections and IUCN data
allThreats <- read_csv(file.path("All_iucn_and_threats.csv")) %>%
  mutate(GroupForPlot = case_when(GroupForPlot == "\nVascular\nplants" ~ "Vascular\nplants",
                                  GroupForPlot == "Terrestrial\n& freshwater\nvertebrates" ~ "Terrestrial &\nfreshwater vertebrates",
                                  GroupForPlot == "Terrestrial\n& freshwater\ninvertebrates" ~ "Terrestrial &\nfreshwater invertebrates",
                                  GroupForPlot == "\n\nMarine" ~ "\nMarine")) %>%
  mutate(GroupForPlot = factor(GroupForPlot, levels = c("Vascular\nplants",
                                                        "Terrestrial &\nfreshwater vertebrates",
                                                        "Terrestrial &\nfreshwater invertebrates",
                                                        "\nMarine")))

### Remove species labelled as extinct
allThreats <- allThreats %>%
  filter(redlistCategory != "EX")

### Taxon group information
taxonInfo  <- read_csv(file.path(projDir, "Global_totals.csv")) %>%
  select(Taxon, Group, Taxon_Higher, Type, Name, Colour) %>%
  mutate(Group = gsub(" NA", "", paste(Group, Taxon_Higher, sep = " "))) %>%
  mutate(GroupForPlot = gsub(" ", "\n", Group, fixed = TRUE)) %>%
  mutate(GroupForPlot = gsub("and\n", "& ", GroupForPlot, fixed = TRUE)) %>%
  mutate(TaxonForPlot = case_when(Name == "phasmids" ~ "Stick and leaf insects",
                                  Name == "sharks"   ~ "Sharks and rays",
                                  .default = Taxon))

### Bioregion information
regionInfo <- read_csv(file.path(projDir, "Bioregion_info.csv")) %>%
  mutate(LabelForPlot = gsub(" ", "\n", Label, fixed = TRUE))

### BTAS map for plotting
base_map        <- readRDS(file.path(mapDir, "gadm_base_map.rds"))
terr_bioregions <- readRDS(file.path(mapDir, "terr_bioregions.rds"))
mare_bioregions <- readRDS(file.path(mapDir, "mare_bioregions.rds")) %>%
  st_crop(st_bbox(base_map))

### Get dividing lines between subregions with land-sea border removed
base_lines <- st_union(st_cast(select(base_map, geom), "MULTILINESTRING"))
mare_lines <- st_union(st_cast(select(mare_bioregions, geom), "MULTILINESTRING"))
terr_lines <- st_union(st_cast(select(terr_bioregions, geom), "MULTILINESTRING"))

base_map_buffer <- st_buffer(base_lines, 10000)
terr_divides <- st_difference(terr_lines, base_map_buffer)
mare_divides <- st_difference(mare_lines, base_map_buffer)

### Colour scheme for Red List categories
colsRedList <- c(#"EX"             = "black",
                 "CR"             = "#A91F13", # rgb(169,  31, 19, maxColorValue = 255), 
                 "EN"             = "orange2",      # rgb(238, 154,  0, maxColorValue = 255),
                 "VU"             = "#F9DA14", # rgb(249, 218, 20, maxColorValue = 255),
                 "Non-threatened" = "grey80",
                 "DD"             = "grey60",
                 "Not evaluated"  = "grey40")

####################################################################################################
### Assign threat types
####################################################################################################

### Convert threats info (which is list of years) to 1s and 0s
threatLevels <- allThreats %>%
  mutate(across(`1_1`:`12_1`, ~ case_when(!is.na(.) ~ 1,
                                           is.na(.) ~ 0)))

### Generate new column for if a species has no listed threat and if it is not assessed
threatLevels <- threatLevels %>%
  mutate(Not_evaluated = case_when(redlistCategory == "Not evaluated" ~ 1,
                                   redlistCategory != "Not evaluated" ~ 0)) %>%
  mutate(No_threats = rowSums(.[,  which(names(allThreats) == "1_1"):which(names(allThreats) == "12_1")])) %>%
  mutate(No_threat_listed = case_when(No_threats == 0 & Not_evaluated == 0 ~ 1,
                                      No_threats >= 1 | Not_evaluated == 1 ~ 0)) %>%
  select(-No_threats)

### Assign threat groups according to Hogue & Breon (2022)
subregionThreatLevels <- threatLevels %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio, names_to = "Subregion", values_to = "Presence") %>%
  relocate(Subregion) %>%
  filter(Presence == 1) %>%
  pivot_longer(cols = "1_1":"No_threat_listed",
               names_to = "ThreatLevel",
               values_to = "Threat_info") %>%
  filter(Threat_info == 1) %>%
  select(-Threat_info) %>%
  mutate(ThreatGroup = case_match(ThreatLevel,
                                  c("1", "1_1", "1_2", "1_3", "2_1_1", "2_1_2", "2_1_3", "2_1_4", "2_2_1",
                                    "2_2_2", "2_2_3", "2_3_1", "2_3_2", "2_3_3", "2_3_4", "2_4_1",  
                                    "2_4_2", "2_4_3", "3_1", "3_2", "3_3", "4_1", "4_2", "4_3", "5_3_3",
                                    "5_3_4", "5_3_5", "6_1", "6_2", "6_3", "7_1_1", "7_1_2", "7_1_3",
                                    "7_2_1", "7_2_10", "7_2_11", "7_2_2", "7_2_3", "7_2_4", "7_2_5",
                                    "7_2_6", "7_2_7", "7_2_8", "7_2_9", "7_3") ~ "Habitat destruction",
                                  c("5_1_1", "5_1_2", "5_1_3", "5_1_4", "5_2_1", "5_2_2", "5_2_3", "5_2_4",
                                    "5_3_1", "5_3_2", "5_4_1", "5_4_2", "5_4_3", "5_4_4", "5_4_5",
                                    "5_4_6")                                   ~ "Overexploitation",
                                  c("8_1_1", "8_1_2", "8_2", "8_2_1", "8_2_2", "8_3", "8_4_1", "8_4_2",
                                    "8_5_1", "8_5_2", "8_6")                   ~ "Invasive species, genes & diseases",
                                  c("9_1_1", "9_1_2", "9_1_3", "9_2_1", "9_2_2", "9_2_3", "9_3_1",
                                    "9_3_2", "9_3_3", "9_3_4", "9_4", "9_5_1", "9_5_2", "9_5_3", "9_5_4",
                                    "9_6_1", "9_6_2", "9_6_3", "9_6_4")        ~ "Pollution",
                                  c("11_1", "11_2", "11_3", "11_4", "11_5")    ~ "Climate change & weather",
                                  c("No_threat_listed")                        ~ "No_threat_listed",
                                  c("Not_evaluated")                            ~ "Not_evaluated")) %>%
  relocate(c(ThreatLevel, ThreatGroup), .after = AsiaEndemic)

### Not all threats are assigned to threat groups in Hogue & Breon
filter(subregionThreatLevels, !is.na(ThreatLevel) & is.na(ThreatGroup)) %>%
  select(ThreatLevel, ThreatGroup) %>%
  distinct()

### Add to new group
subregionThreatLevels <- subregionThreatLevels %>%
  mutate(ThreatGroup = case_when(ThreatLevel %in% c("4_4",
                                                    "10_1", "10_2", "10_3",
                                                    "12_1") ~ "Natural disasters & others",
                                 .default = ThreatGroup)) %>%
  mutate(ThreatGroup = factor(ThreatGroup,
                              levels = rev(c("Habitat destruction",
                                             "Overexploitation",
                                             "Invasive species, genes & diseases",
                                             "Pollution",
                                             "Climate change & weather",
                                             "Natural disasters & others",
                                             "No_threat_listed",
                                             "Not_evaluated")),
                              labels = rev(c("Habitat destruction",
                                             "Overexploitation",
                                             "Invasive species,\ngenes & diseases",
                                             "Pollution",
                                             "Climate change\n& weather",
                                             "Natural disasters\n& others",
                                             "No threat listed",
                                             "Not evaluated")))) %>%
  mutate(GroupForPlot = factor(GroupForPlot, levels = c("Vascular\nplants",
                                                        "Terrestrial &\nfreshwater vertebrates",
                                                        "Terrestrial &\nfreshwater invertebrates",
                                                        "\nMarine")))

### Summarise at taxon-level
taxonThreatLevels <- reframe(subregionThreatLevels,
                             Species      = unique(Species),
                             TaxonForPlot = unique(TaxonForPlot),
                             GroupForPlot = unique(GroupForPlot),
                             .by = c("UniqueName", "ThreatGroup")) %>%
  mutate(ThreatGroup = factor(ThreatGroup,
                              levels = rev(c("Habitat destruction",
                                             "Overexploitation",
                                             "Invasive species,\ngenes & diseases",
                                             "Pollution",
                                             "Climate change\n& weather",
                                             "Natural disasters\n& others",
                                             "No threat listed",
                                             "Not evaluated")),
                              labels = rev(c("Habitat destruction",
                                             "Overexploitation",
                                             "Invasive species,\ngenes & diseases",
                                             "Pollution",
                                             "Climate change\n& weather",
                                             "Natural disasters\n& others",
                                             "No threat listed",
                                             "Not evaluated")))) %>%
  mutate(GroupForPlot = factor(GroupForPlot, levels = c("Vascular\nplants",
                                                        "Terrestrial &\nfreshwater vertebrates",
                                                        "Terrestrial &\nfreshwater invertebrates",
                                                        "\nMarine")))

####################################################################################################
### Plot: Prop threatened for tropical Asia - by group
####################################################################################################

### Get regional richness for pie size
regRich <- allThreats %>%
  reframe(RichRegion = n(), .by = GroupForPlot)

### Summarise each taxon
propThreatAll <- allThreats %>%
  select(GroupForPlot, Species, redlistCategory) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, redlistCategory)) %>%
  mutate(redlistCategory = factor(redlistCategory,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(regRich, by = "GroupForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(LogRichRegion = log(RichRegion)) %>%
  mutate(Subset = "All species")

### Summary of number of threatened species as proportion of assessed species
propThreatAll %>% filter(redlistCategory != "Not evaluated") %>% 
  reframe(nEval   = sum(NoThreat),
          nThreat = sum(NoThreat[redlistCategory %in% c("CR", "EN", "VU")]),
          prop    = sum(NoThreat[redlistCategory %in% c("CR", "EN", "VU")]) / sum(NoThreat),
          CR      = NoThreat[redlistCategory == "CR"],
          EN      = NoThreat[redlistCategory == "EN"],
          VU      = NoThreat[redlistCategory == "VU"],
          .by = GroupForPlot)

### Repeat but with non-assessed species removed
### Get regional richness for pie size
regRich <- allThreats %>%
  filter(redlistCategory != "Not evaluated") %>%
  reframe(RichRegion = n(), .by = GroupForPlot)

### Summarise each taxon
propThreatEval <- allThreats %>%
  filter(redlistCategory != "Not evaluated") %>%
  select(GroupForPlot, Species, redlistCategory) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, redlistCategory)) %>%
  mutate(redlistCategory = factor(redlistCategory,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(regRich, by = "GroupForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(LogRichRegion = log(RichRegion)) %>%
  mutate(Subset = "Evaluated\nspecies only")

### Combine
propThreat <- bind_rows(propThreatAll, propThreatEval) %>%
  mutate(Group = paste(GroupForPlot, Subset))

### Loop through and plot
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

### Create final plot
pIucnCatsGroup <-
  propThreat_group_plot_list[[1]] + propThreat_group_plot_list[[2]] + 
  propThreat_group_plot_list[[3]] + propThreat_group_plot_list[[4]] +
  propThreat_group_plot_list[[5]] + propThreat_group_plot_list[[6]] + 
  propThreat_group_plot_list[[7]] + propThreat_group_plot_list[[8]] +
  plot_layout(guides = "collect", tag_level = "keep", nrow = 2, ncol = 4, axis_titles = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_blank())
pIucnCatsGroup

ggsave(file.path(figDir, "SI", "Prop_threatened_by_group.png"), pIucnCatsGroup,
       width = 180, height = 110, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "SI", "Prop_threatened_by_group.pdf"), pIucnCatsGroup,
       width = 180, height = 110, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### Plot: Prop threatened by subregion - map version - All taxa
####################################################################################################

### Get subregions richness for standardisation
subregRich <- allThreats %>%
  reframe(across(Indian_Subcontinent:South_Kuroshio, \(x) sum(x, na.rm = TRUE)), .by = GroupForPlot)  %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio,
               names_to = "Subregion",
               values_to = "RichSubregion")

### Get subregion centroids
coords <- bind_rows(st_drop_geometry(terr_bioregions),
                    rename(st_drop_geometry(mare_bioregions), Bioregion = PROVINCE)) %>%
  mutate(Bioregion = gsub(" ", "_", Bioregion)) %>%
  select(Bioregion, X, Y) %>%
  mutate(X = case_when(Bioregion == "Sulawesi" ~ 13490000,
                       Bioregion == "Andaman"  ~ 10400000,
                       Bioregion == "Bay_of_Bengal"                ~  9600000,
                       Bioregion == "Central_Indian_Ocean_Islands" ~  8150000,
                       Bioregion == "Eastern_Coral_Triangle"       ~ 16700000,
                       Bioregion == "Sahul_Shelf"                  ~ 15100000,
                       Bioregion == "South_China_Sea"              ~ 12800000,
                       Bioregion == "West_and_South_Indian_Shelf"  ~  9000000,
                       .default = X)) %>%
  mutate(Y = case_when(Bioregion == "Sumatra"  ~ -370000,
                       Bioregion == "Sulawesi" ~ -180000,
                       Bioregion == "Andaman"  ~ 1100000,
                       Bioregion == "Eastern_Coral_Triangle"      ~ -600000,
                       Bioregion == "South_China_Sea"             ~ 1750000,
                       Bioregion == "West_and_South_Indian_Shelf" ~  760000,
                       .default = Y))

### Summarise across subregions
subregionIucn <- allThreats %>%
  pivot_longer(cols = Indian_Subcontinent:South_Kuroshio, names_to = "Subregion", values_to = "Presence") %>%
  relocate(Subregion) %>%
  filter(Presence == 1) %>%
  select(GroupForPlot, Subregion, redlistCategory) %>%
  reframe(NoThreat = n(),
          .by = c(GroupForPlot, Subregion, redlistCategory)) %>%
  mutate(redlistCategory = factor(redlistCategory,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR"))) %>%
  left_join(subregRich, by = c("GroupForPlot", "Subregion")) %>%
  mutate(propThreat = NoThreat / RichSubregion) %>%
  mutate(LogRichSubregion = log(RichSubregion)) %>%
  left_join(select(regionInfo, Bioregion), join_by("Subregion" == "Bioregion")) %>%
  left_join(coords, join_by("Subregion" == "Bioregion")) %>%
  pivot_wider(id_cols = c(GroupForPlot, Subregion, LogRichSubregion, X, Y),
              names_from = redlistCategory,
              values_from = propThreat)
subregionIucn[is.na(subregionIucn)] <- 0

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

### Add in subplots
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
  inset_element(pBox,       left = 0.71, right = 0.98, bottom = 0.62, top = 0.94) +
  inset_element(pSubPlants, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubVerts <- 
  propThreat_group_plot_list[["All_verts"]] + propThreat_group_plot_list[["Eval_verts"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[2]] <- region_iucn_plot_list[[2]] + 
  inset_element(pBox,      left = 0.71, right = 0.98, bottom = 0.62, top = 0.94) +
  inset_element(pSubVerts, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubInverts <- 
  propThreat_group_plot_list[["All_inverts"]] + propThreat_group_plot_list[["Eval_inverts"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[3]] <- region_iucn_plot_list[[3]] + 
  inset_element(pBox,        left = 0.71, right = 0.98, bottom = 0.62, top = 0.94) +
  inset_element(pSubInverts, left = 0.71, right = 0.98, bottom = 0.59, top = 0.98)

pSubMarine <-
  propThreat_group_plot_list[["All_marine"]] + propThreat_group_plot_list[["Eval_marine"]] & 
  theme(plot.title      = element_text(family = "arial", colour = "black", size = 5, hjust = 0.5), 
        strip.text.x    = element_blank(),
        strip.text.y    = element_blank(),
        legend.position = "none")
region_iucn_plot_list[[4]] <- region_iucn_plot_list[[4]] +
  inset_element(pBoxMarine, left = 0.72, right = 0.99, bottom = 0.48, top = 0.79) +
  inset_element(pSubMarine, left = 0.72, right = 0.99, bottom = 0.45,  top = 0.84)

pIucnCatsSubregionMap <-
  (region_iucn_plot_list[[1]] +  region_iucn_plot_list[[2]]) /
  (region_iucn_plot_list[[3]] +  region_iucn_plot_list[[4]]) +
  plot_layout(tag_level = "new", guides = "collect") &
  theme(legend.position = "bottom")
pIucnCatsSubregionMap

ggsave(pIucnCatsSubregionMap,
       filename = file.path(figDir, "map_prop_threatened.png"),
       width = 180, height = 110, units = "mm", dpi = 600, bg = "white")

ggsave(pIucnCatsSubregionMap,
       filename = file.path(figDir, "map_prop_threatened.pdf"),
       width = 180, height = 110, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### Plot: Heat map plots
####################################################################################################

### Summary of number of threatened species as proportion of assessed species
taxonThreatLevels %>% 
  filter(!ThreatGroup %in% c("No threat listed", "Not evaluated")) %>% 
  reframe(nEval    = n(),
          propEval = n() / length(unique(.$UniqueName)),
          .by = c(ThreatGroup))

### Summary of number of threatened species as proportion of assessed species by taxonomic group
taxonThreatLevels %>%
  filter(!ThreatGroup %in% c("No threat listed", "Not evaluated")) %>%
  group_by(GroupForPlot) %>%
  mutate(nEval = length(unique(UniqueName))) %>%
  ungroup() %>%
  reframe(nEval      = unique(nEval),
          nThreat    = n(),
          propThreat = n() / unique(nEval),
          .by = c(ThreatGroup, GroupForPlot)) %>%
  pivot_wider(id_cols = GroupForPlot,
              names_from = ThreatGroup,
              values_from = propThreat)
  
### Get regional richness for pie size
regRich <- taxonThreatLevels %>%
  filter(!ThreatGroup %in% c("No threat listed", "Not evaluated")) %>%
  reframe(RichRegion = length(unique(Species)), .by = TaxonForPlot)

### Calculate as proportion of taxon richness
propType <- taxonThreatLevels %>%
  filter(!ThreatGroup  %in% c("No threat listed", "Not evaluated")) %>%
  select(TaxonForPlot, Species, ThreatGroup) %>%
  distinct() %>%
  reframe(NoThreat = length(unique(Species)),
          .by = c(TaxonForPlot, ThreatGroup))

### Add in combos with no species (set prop to 0)
missingEntries <- expand.grid(TaxonForPlot = unique(subregionThreatLevels$TaxonForPlot),
                              ThreatGroup  = unique(subregionThreatLevels$ThreatGroup),
                              NoThreat     = NA) %>%
  filter(!ThreatGroup  %in% c("No threat listed", "Not evaluated")) %>%
  filter(!paste(TaxonForPlot, ThreatGroup) %in% paste(propType$TaxonForPlot, propType$ThreatGroup))

propType <- propType %>%
  bind_rows(missingEntries) %>%
  full_join(regRich, by = "TaxonForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(propThreat = case_when(is.na(propThreat) ~ 0, .default = propThreat)) %>%
  mutate(propThreat = case_when(TaxonForPlot %in% c("Ants", "Bees", "Caddisflies", "Centipedes",
                                                    "Mirid bugs", "Sponges") ~ NA,
                                .default = propThreat)) %>%
  full_join(taxonInfo, by = "TaxonForPlot") %>%
  mutate(TaxonForPlot = factor(TaxonForPlot,
                               levels = taxonInfo$Taxon))

pThreatsHeat <- ggplot(propType, aes(TaxonForPlot, ThreatGroup)) +
  theme(panel.background = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", family = "sans", size = 6),
        axis.text.y = element_text(colour = "black", family = "sans", size = 6),
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"),
        legend.position       = "bottom",
        legend.justification  = "center",
        legend.key.width      = unit(10, "mm"),
        legend.key.height     = unit(3, "mm"),
        legend.frame          = element_rect(colour = "black", linewidth = 0.25),
        legend.ticks          = element_line(colour = "black", linewidth = 0.25),
        legend.ticks.length   = unit(c(-1, 0), "mm"),
        legend.title.position = "left",
        legend.box.margin     = margin(t = -4, r = 0, b = 0, l = 0, unit = "mm"),
        legend.title          = element_text(family = "sans", face = "plain", colour = "black", size = 6, vjust = 1, hjust = 0.5,
                                             margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm")),
        legend.text           = element_text(family = "sans", colour = "black", size = 5)
  ) +
  labs(x = NULL, y = NULL) +# y = "Threat") +
  geom_tile(aes(fill = propThreat), col = "grey70") +
  scale_fill_distiller(palette = "OrRd", limits = c(0, 1), direction = 1,
                       name = "Prop. of species with\nthreat information") +
  coord_fixed(expand = FALSE, clip = "off") +
  geom_rect(xmin = 0.5,  xmax =  4.5, ymin = -Inf, ymax = 8, colour = "black", fill = NA, linewidth = 0.25) +
  geom_rect(xmin = 4.5,  xmax =  9.5, ymin = -Inf, ymax = 8, colour = "black", fill = NA, linewidth = 0.25) +
  geom_rect(xmin = 9.5,  xmax = 20.5, ymin = -Inf, ymax = 8, colour = "black", fill = NA, linewidth = 0.25) +
  geom_rect(xmin = 20.5, xmax = 24.5, ymin = -Inf, ymax = 8, colour = "black", fill = NA, linewidth = 0.25) +
  geom_rect(xmin = 0.5,  xmax = 24.5, ymin = -Inf, ymax = Inf, colour = "black", fill = NA, linewidth = 0.25) +
  geom_text(x = 2.5,  y = 7.2, label = "Vascular\nplants",
            family = "sans", hjust = 0.5, vjust = 0.5, fontface = "plain", size.unit = "pt", size = 8, colour = "black") +
  geom_text(x = 7,    y = 7.2, label = "Terrestrial &\nfreshwater vertebrates",
            family = "sans", hjust = 0.5, vjust = 0.5, fontface = "plain", size.unit = "pt", size = 8, colour = "black") +
  geom_text(x = 15,   y = 7.2, label = "Terrestrial &\nfreshwater invertebrates",
            family = "sans", hjust = 0.5, vjust = 0.5, fontface = "plain", size.unit = "pt", size = 8, colour = "black") +
  geom_text(x = 22.5, y = 7.2, label = "Marine",
            family = "sans", hjust = 0.5, vjust = 0.5, fontface = "plain", size.unit = "pt", size = 8, colour = "black")
pThreatsHeat

ggsave(pThreatsHeat,
       filename = file.path(figDir, "Threats_heat_map.png"),
       width = 180, height = 85, units = "mm", dpi = 600, bg = "white")

ggsave(pThreatsHeat,
       filename = file.path(figDir, "Threats_heat_map.pdf"),
       width = 180, height = 85, units = "mm", dpi = 600, bg = "white")

####################################################################################################
### Plot: Year published - by group
####################################################################################################

for(i in 1:length(unique(allThreats$GroupForPlot))) {
  ### Subset taxon
  taxon <- unique(allThreats$GroupForPlot)[i]
  dd <- droplevels(filter(allThreats, GroupForPlot == taxon)) %>%
    filter(!is.na(yearPublished))
  
  ### If there are some IUCN data for taxon
  if(nrow(dd) > 0) {
    p <- ggplot(dd, aes(yearPublished)) +
      theme(panel.background = element_rect(fill = NA, colour = "black"),
            title = element_text(family = "sans", size = 8, colour = "black", hjust = 0.5),
            axis.text.x = element_text(family = "sans", size = 6, colour = "black"),#, angle = 90, vjust = 0.5),
            axis.text.y = element_text(family = "sans", size = 6, colour = "black"),
            axis.title  = element_text(family = "sans", size = 8, colour = "black"),
            axis.ticks  = element_line(colour = "black", linewidth = 0.25),
            axis.line   = element_line(colour = "black", linewidth = 0.25),
            strip.text  = element_text(family = "sans", size = 8, colour = "black"),
            strip.background = element_blank(),# element_rect(fill = NA, colour = "black")) +
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey50", linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = "grey50", linewidth = 0.1)) +
      facet_wrap(GroupForPlot ~ .) +#, switch = "y") +
      lims(x = c(1995, 2025)) +
      labs(x = "Assessment year",
           y = "No. of species") +
      scale_x_continuous(breaks = function(x) unique(floor(pretty(x, 3))),
                         expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(breaks = function(x) unique(floor(pretty(x, 4))),
                         expand = expansion(mult = c(0, 0.05)),
                         labels = label_number(accuracy = 1, big.mark = "")) +
      geom_histogram(breaks = seq(1995.5, 2025.5, 1))
  }
  assign(paste0("p", i), p)
}
pAssessmentYearGroup <- p1 + p2 + p3 + p4 +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect_x",
              nrow = 1, ncol = 4)
# pAssessmentYearGroup

ggsave(file.path(figDir, "Assessment_year.png"), pAssessmentYearGroup,
       width = 180, height = 50, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "Assessment_year.pdf"), pAssessmentYearGroup,
       width = 180, height = 50, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### Figure 6: Put all three parts together
####################################################################################################

threatsCombined <-
  (pIucnCatsSubregionMap & theme(legend.position = "bottom")) +
  pThreatsHeat +
  pAssessmentYearGroup + 
  plot_layout(guides = "keep", tag_level = "keep", 
              axis_titles = "keep", heights = c(0.3, 0.3, 0.2, 0.15)) +
  plot_annotation(title = c("A", "B", "C"))

  
ggsave(file.path(figDir, "Threats_combined.png"), threatsCombined,
       width = 180, height = 250, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "Threats_combined.pdf"), threatsCombined,
       width = 180, height = 230, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### Plot: Prop threatened for tropical Asia - each taxon separately
####################################################################################################

### Get regional richness for pie size
regRich <- allThreats %>%
  reframe(RichRegion = n(), .by = TaxonForPlot)

### Summarise each taxon
propThreat <- allThreats %>%
  select(TaxonForPlot, GroupForPlot, Species, redlistCategory) %>%
  reframe(NoThreat = n(),
          .by = c(TaxonForPlot, GroupForPlot, redlistCategory)) %>%
  mutate(redlistCategory = factor(redlistCategory,
                                  levels = c("Not evaluated", "DD", "Non-threatened", "VU","EN", "CR", "EX"))) %>%
  left_join(regRich, by = "TaxonForPlot") %>%
  mutate(propThreat = NoThreat / RichRegion) %>%
  mutate(LogRichRegion = log(RichRegion))

for(i in 1:length(unique(propThreat$TaxonForPlot))) {
  ### Subset taxon
  taxon <- unique(allThreats$TaxonForPlot)[i]
  dd <- filter(propThreat, TaxonForPlot == taxon)
  
  ### Generate plot
  p <- ggplot(dd,
              aes(x = LogRichRegion / 2,
                  y = NoThreat,
                  fill = redlistCategory,
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
  p21 + p22 + p23 + p24 +
  plot_layout(guides = "collect", tag_level = "keep", nrow = 5, ncol = 6, axis_titles = "collect") &
  theme(legend.position = "bottom")
pIucnCatsAll

ggsave(file.path(figDir, "SI", "Prop_threatened.png"), pIucnCatsAll,
       width = 180, height = 175, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "SI", "Prop_threatened.pdf"), pIucnCatsAll,
       width = 180, height = 175, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

####################################################################################################
### Plot: Year published - each taxon separately
####################################################################################################

for(i in 1:length(unique(allThreats$TaxonForPlot))) {
  ### Subset taxon
  taxon <- unique(allThreats$TaxonForPlot)[i]
  dd <- droplevels(filter(allThreats, TaxonForPlot == taxon)) %>%
    filter(!is.na(yearPublished))
  
  ### If there are some IUCN data for taxon
  if(nrow(dd) > 0) {
    p <- ggplot(dd, aes(yearPublished)) +
      theme(panel.background = element_rect(fill = NA, colour = "black"),
            title = element_text(family = "sans", size = 8, colour = "black", hjust = 0.5),
            axis.text.x = element_text(family = "sans", size = 6, colour = "black"),#, angle = 90, vjust = 0.5),
            axis.text.y = element_text(family = "sans", size = 6, colour = "black"),
            axis.title  = element_text(family = "sans", size = 8, colour = "black"),
            axis.ticks  = element_line(colour = "black", linewidth = 0.25),
            axis.line   = element_line(colour = "black", linewidth = 0.25),
            strip.text  = element_text(family = "sans", size = 6, colour = "black"),
            strip.background = element_blank(),# element_rect(fill = NA, colour = "black")) +
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey50", linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = "grey50", linewidth = 0.1)) +
      facet_wrap(TaxonForPlot ~ ., switch = "y") +
      lims(x = c(1996, 2025)) +
      scale_x_continuous(breaks = function(x) unique(floor(pretty(x, 3))),
                         expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(breaks = function(x) unique(floor(pretty(x, 4))), position = "right",
                         expand = expansion(mult = c(0, 0.05)),
                         labels = label_number(accuracy = 1, big.mark = "")) +
      geom_histogram(breaks = seq(1995.5, 2025.5, 1))
  }
  
  ### If there is no IUCN data generate blank plot
  if(nrow(dd) == 0) {
    dd <- data.frame(TaxonForPlot = taxon)
    p <- ggplot(dd) +
      theme(panel.background = element_rect(fill = NA, colour = "black"),
            title = element_text(family = "sans", size = 8, colour = "black", hjust = 0.5),
            axis.text.x = element_text(family = "sans", size = 6, colour = "black"),#, angle = 90, vjust = 0.5),
            axis.text.y = element_text(family = "sans", size = 6, colour = "black"),
            axis.title  = element_text(family = "sans", size = 8, colour = "black"),
            axis.ticks  = element_line(colour = "black", linewidth = 0.25),
            axis.line   = element_line(colour = "black", linewidth = 0.25),
            strip.text  = element_text(family = "sans", size = 6, colour = "black"),
            strip.background = element_blank(),# element_rect(fill = NA, colour = "black")) +
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey50", linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = "grey50", linewidth = 0.1)) +
      facet_wrap(TaxonForPlot ~ ., switch = "y") +
      scale_y_continuous(position = "right") +
      lims(x = c(1995.5, 2025.5))
  }
  
  ### Add axis labels to specific plots
  if(i %in% c(1, 5, 10, 16, 21)) {
    p <- p +
      facet_grid(TaxonForPlot ~ GroupForPlot, switch = "y")
  }
  if(i %in% c(15)) {#c(4, 9, 15, 20, 24)) {
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
  p4 +            p8 +            p13 +            p19 +           p24 +
  plot_spacer() + p9 +            p14 +            p20 +           plot_spacer() +
  plot_spacer() + plot_spacer() + p15 +            plot_spacer() + plot_spacer() +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect_x", tag_level = "keep",
              nrow = 6, ncol = 5)
# pAssessmentYear

ggsave(file.path(figDir, "SI", "Assessment_year_by_taxon.png"), pAssessmentYear,
       width = 180, height = 155, units = "mm", dpi = 600, bg = "white")

ggsave(file.path(figDir, "SI", "Assessment_year_by_taxon.pdf"), pAssessmentYear,
       width = 220, height = 155, units = "mm", dpi = 600, bg = "white", device = cairo_pdf)

