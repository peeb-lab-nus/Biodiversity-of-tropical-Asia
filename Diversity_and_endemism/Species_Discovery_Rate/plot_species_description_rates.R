# GENERATE SPECIES DESCRIPTIONS CURVES ============

rm(list = ls())

## PACKAGES ====================
library(plyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)

## DIRECTORIES ====================

# Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
proj.dir <- "Diversity_and_endemism"

# You shouldn't need to adjust these folders
data.dir <- file.path(proj.dir, "Endemicity")
fun.dir  <- file.path(proj.dir, "Analysis_functions")
res.dir  <- file.path(proj.dir, "Species_Discovery_Rate")
fig.dir  <- file.path(proj.dir, "Species_Discovery_Rate", "Figures")
if(!dir.exists(fig.dir)) { dir.create(fig.dir, recursive = TRUE) }

## IMPORT DATA ====================

### Load function for geernating discovery curves
source(file.path(fun.dir, "Discovery_rates", "generate_description_curves.R"))

### Spreadsheet with taxon info
taxonInfo <- read.csv(file.path(proj.dir, "Global_totals.csv")) %>%
  arrange(desc(Group), desc(Taxon_Higher), Taxon)

### Read in endemicity data (out-putted from calc_endemicity.R)
taxa_full_list <- read.csv(file.path(data.dir, "taxon_full_list.csv"))
nrow(taxa_full_list)

sum(is.na(taxa_full_list$Year)) # 3 (2 plants and 1 millipede)
sum(is.infinite(taxa_full_list$Year)) # 0
subset(taxa_full_list, is.na(Year))

taxa_full_list_clean <- subset(taxa_full_list, !(is.na(Year)| is.infinite(Year)))
length(unique(taxa_full_list_clean$Taxon)) # 24 taxa
range(taxa_full_list_clean$Year) # 1753 to 2024

table(taxa_full_list[c("Taxon", "Region")]) # Some species with unknown region

### Read in bioregion info which has info on plot order and colours etc
bioregion_info <- read.csv(file.path(proj.dir, "Bioregion_info.csv")) %>%
  filter(!is.na(Order)) %>%
  mutate(Bioregion = gsub(" ", "_", Bioregion)) %>%
  mutate(Label = case_when(Label == "Malay peninsula"              ~ "Malay\nPeninsula",
                           Label == "Lesser Sunda Islands"         ~ "Lesser Sunda\nIslands",
                           Label == "Central Indian Ocean Islands" ~ "Central Indian\nOcean Islands",
                           Label == "West and South Indian Shelf"  ~ "West & South\nIndian Shelf",
                           Label == "Bay of Bengal"                ~ "Bay of\nBengal",
                           Label == "Java Transitional"            ~ "Java\nTransitional",
                           Label == "South China Sea"              ~ "South\nChina Sea",
                           Label == "South Kuroshio"               ~ "South\nKuroshio",
                           Label == "Western Coral Triangle"       ~ "Western\nCoral Triangle",
                           Label == "Eastern Coral Triangle"       ~ "Eastern\nCoral Triangle",
                           .default = Label))

# ==================================================================================================
## GENERATE SUMMARY STATISTICS
# ==================================================================================================

# Define regions
terr_regions         <- c(bioregion_info$Bioregion[bioregion_info$Realm == "Terrestrial"], "Multiple regions")
marine_regions       <- c(bioregion_info$Bioregion[bioregion_info$Realm == "Marine"], "Multiple regions")
terr_region_labels   <- c(bioregion_info$Label[bioregion_info$Realm == "Terrestrial"], "Multiple regions")
marine_region_labels <- c(bioregion_info$Label[bioregion_info$Realm == "Marine"], "Multiple regions")
region_colors        <- c(bioregion_info$Colour[bioregion_info$Realm == "Terrestrial"], "grey80")

combined_regions       <- c(bioregion_info$Bioregion, "Multiple regions")
combined_region_labels <- c(bioregion_info$Label, "Multiple regions")
combined_region_colors <- c(bioregion_info$Colour, "grey80")

# Define taxa
taxonInfo <- taxonInfo %>%
  mutate(Order = 1:n()) %>%
  mutate(FacetGroup = case_when(Group == "Vascular plants" ~ "Vascular plants",
                                Group == "Terrestrial and freshwater" &
                                  Taxon_Higher == "vertebrates" ~ "Terrestrial & freshwater\nvertebrates",
                                Group == "Terrestrial and freshwater" & 
                                  Taxon_Higher == "invertebrates" ~ "Terrestrial & freshwater\ninvertebrates",
                                Group == "Marine" ~ "Marine"))

taxon_colors <- taxonInfo$Colour
taxon_order  <- taxonInfo$Taxon
facet_groups <- unique(unique(taxonInfo$FacetGroup))
facet_group_abbrev <- c("plants",
                        "terr_fw_vertebrates",
                        "terr_fw_invertebrates",
                        "marine")

# ==================================================================================================
## GENERATE SUMMARY STATISTICS
# ==================================================================================================

# Get range for each taxon
taxonomic_summary_statistics <- ddply(.data = taxa_full_list_clean,
                                      .variables = .(Taxon),
                                      .fun = summarise,
                                      minYear = min(Year),
                                      maxYear = max(Year), 
                                      S_2011_2020 = length(unique(Species[Year <= 2020 & Year > 2010])),
                                      S_2001_2010 = length(unique(Species[Year <= 2010 & Year > 2000])),
                                      S_1991_2000 = length(unique(Species[Year <= 2000 & Year > 1990])))
write.csv(taxonomic_summary_statistics,
          file.path(res.dir, "taxonomic_summary_statistics.csv"),
          row.names = FALSE)

# ==================================================================================================
## SPECIES ACCUMULATION FOR WHOLE-REGION
# ==================================================================================================

### Extract discovery curves by taxon
taxa_description_curves_df <- ddply(.data = taxa_full_list_clean,
                                    .variables = .(Taxon),
                                    .fun = function(x){
                                      x %>%
                                        .[,c("Species", "Year")] %>%
                                        # Take unique so it is at region level
                                        unique() %>%
                                        .[["Year"]] %>%
                                        generate_description_curves()
                                    })

### Scale
taxa_description_curves_scaled_df <- ddply(.data = taxa_description_curves_df,
                                           .variables = .(Taxon),
                                           .fun = function(x){
                                             x$S_scaled <- x$S / max(x$S)
                                             return(x)
                                           })

taxa_description_curves_scaled_df <- left_join(taxa_description_curves_scaled_df, taxonInfo)
taxa_description_curves_scaled_df$Taxon      <- factor(taxa_description_curves_scaled_df$Taxon,
                                                       levels = taxon_order)
taxa_description_curves_scaled_df$FacetGroup <- factor(taxa_description_curves_scaled_df$FacetGroup,
                                                       levels = facet_groups)

### Generate plots
discovery_accumulation_alt_plot <- ggplot() +
  geom_path(aes(y = S_scaled,
                x = t,
                group = Taxon,
                linetype = Taxon,
                colour = Taxon
  ),
  data = taxa_description_curves_scaled_df,
  linewidth = 0.5) +
  facet_wrap(~ FacetGroup, nrow = 1) +
  labs(x = "Year", y = "Proportion of species\ndescribed") +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(1750, 1800, 1850, 1900, 1950, 2000),
                     limits = c(1749, 2025)) +
  scale_colour_manual(values = taxon_colors) +
  scale_linetype_manual(values = c(rep(1:2, 2)[1:4],    # Plants
                                   rep(1:2, 3)[1:5],    # Verts
                                   rep(1:2, 6)[1:11],   # Inverts
                                   rep(1:2, 2)[1:4])) + # Marine
  theme(panel.background   = element_blank(),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.background    = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.text    = element_text(family = "Sans", size = 5, colour = "black"),
        axis.title   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.ticks   = element_line(colour = "black", linewidth = 0.5),
        legend.title = element_blank(),# legend.title = element_text(family = "Sans", size = 7, colour = "black", hjust = 0.5),
        legend.text  = element_text(family = "Sans", size = 5, colour = "black"),
        legend.title.position = "top",
        legend.key = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())
discovery_accumulation_alt_plot

# Plot legends separately to be added later on
for(i in 1:length(facet_groups)){
  test_plot <- ggplot() +
    geom_path(aes(y = S_scaled,
                  x = t,
                  group = Taxon,
                  linetype = Taxon,
                  colour = Taxon
    ),
    data = subset(taxa_description_curves_scaled_df,
                  FacetGroup == facet_groups[i])) +
    guides(colour   = guide_legend(ncol = ifelse(grepl(facet_groups[i], pattern = "invertebrates"), 2, 1)),
           linetype = guide_legend(ncol = ifelse(grepl(facet_groups[i], pattern = "invertebrates"), 2, 1))) +
    scale_colour_manual(values = subset(taxonInfo, FacetGroup == facet_groups[i])$Colour) +
    scale_linetype_manual(values = c(rep(1:2, 2)[1:4],    # Plants
                                     rep(1:2, 3)[1:5],    # Verts
                                     rep(1:2, 6)[1:11],   # Inverts
                                     rep(1:2, 2)[1:4])) + # Marine
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.text  = element_text(size = 6, colour = "black"),# family = "Sans"),
          legend.key.spacing.y = unit(0.01, "mm"),
          legend.box.margin = margin(c(0, 0, 0, 0)),
          legend.margin = margin(c(0, 0, 0, 0)),
          legend.key.size = unit(3, "mm"),
          legend.justification = "top",
          legend.justification.left = 1,
          legend.position = "right")
  FacetGroupLeg <- get_legend(test_plot)
  assign(paste("FacetGroupLeg", facet_group_abbrev[i], sep = "_"), FacetGroupLeg)
}

# ==================================================================================================
## GENERATE NO OF DESCRIPTIONS PER YEAR PLOTS - BY REGION
# ==================================================================================================

### List of species per region with description year
taxa_region_list <- taxa_full_list_clean %>% 
  select(all_of(c("Taxon","Species", "Region", "Year"))) %>%
  subset(!is.na(Year)) %>%
  subset(!Region %in% "Unknown") %>% # about 201 taxa
  unique()

### Which species are found across more than one region
taxa_region_count <- taxa_region_list %>%
  reframe(nRegion = length(unique(Region)),
          .by = c(Taxon, Species))

taxa_region_list2 <- left_join(taxa_region_list, taxa_region_count)

taxa_region_list2$Region[taxa_region_list2$nRegion > 1] <- "Multiple regions"

taxa_region_list3 <- left_join(taxa_region_list2, taxonInfo) %>%
  mutate(FacetGroup = factor(FacetGroup,
                             levels = c("Vascular plants",
                                        "Terrestrial & freshwater\nvertebrates",
                                        "Terrestrial & freshwater\ninvertebrates",
                                        "Marine"),
                             labels = c("\nVascular plants",
                                        "Terrestrial & freshwater\nvertebrates",
                                        "Terrestrial & freshwater\ninvertebrates",
                                        "\nMarine"))) %>%
  mutate(Region = factor(Region,
                         levels = combined_regions,
                         labels = combined_region_labels))

### Plot
description_time_region_cont_plot <-
  ggplot(data = taxa_region_list3) + 
  geom_bar(aes(x = Year, fill = Region), stat = "count") +
  facet_grid(FacetGroup ~ ., scale = "free_y") +
  scale_fill_manual(values = combined_region_colors) +
  scale_x_continuous(breaks = seq(1750, 2000, 50),
                     limits = c(1750, 2025),
                     expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
  guides(fill = guide_legend(title = NULL)) +
  labs(y = "No of species descriptions per year") +
  theme(panel.grid         = element_blank(),
        panel.background   = element_blank(),
        panel.border       = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "grey95", linewidth = 0.75),
        panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.75),
        axis.ticks         = element_line(colour = "black", linewidth = 0.5),
        axis.text.x        = element_text(family = "Sans", size = 5, colour = "black", angle = 45, hjust = 1),
        axis.text.y        = element_text(family = "Sans", size = 5, colour = "black"),
        axis.title         = element_text(family = "Sans", size = 8, colour = "black"),
        strip.text         = element_text(family = "Sans", size = 7, colour = "black"),
        strip.background   = element_blank(),
        legend.text        = element_text(family = "Sans", size = 6, colour = "black"),
        legend.position    = "bottom")
description_time_region_cont_plot

# mock figures to generate legends
terr_description_legend_plot <- 
  ggplot(data = subset(taxa_region_list3, FacetGroup == "\nVascular plants")) + 
  geom_bar(aes(x = Year, fill = Region), stat = "count") +
  scale_fill_manual(values = c(combined_region_colors[1:11], "grey90")) +
  guides(fill = guide_legend(title = "Terrestrial & freshwater", nrow = 4, byrow = FALSE)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.text  = element_text(family = "Sans", size = 5, colour = "black"),
        legend.title = element_text(family = "Sans", size = 7, colour = "black", hjust = 0.5),
        legend.key.size = unit(4, "mm"),
        legend.position = "right")

marine_description_legend_plot <- 
  ggplot(data = subset(taxa_region_list3, FacetGroup == "\nMarine")) + 
  geom_bar(aes(x = Year, fill = Region), stat = "count") +
  scale_fill_manual(values = combined_region_colors[12:23] ) +
  guides(fill = guide_legend(title = "Marine", nrow = 4, byrow = FALSE)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.text  = element_text(family = "Sans", size = 5, colour = "black"),
        legend.title = element_text(family = "Sans", size = 7, colour = "black", hjust = 0.5),
        legend.key.size   = unit(4, "mm"),
        legend.position   = "right")

# ==================================================================================================
## GENERATE FINAL FIGURE 4
# ==================================================================================================

# Panel A with separated legends for each taxonomic group
discovery_accumulation_alt_plot_leg <-
  plot_grid(discovery_accumulation_alt_plot + theme(legend.position = "none"),
            plot_grid(NULL,
                      FacetGroupLeg_plants,
                      FacetGroupLeg_terr_fw_vertebrates,
                      FacetGroupLeg_terr_fw_invertebrates,
                      FacetGroupLeg_marine,
                      ncol = 5, align = "h", rel_widths = c(0.3, 1, 1, 1, 1)),
            nrow = 2, rel_heights = c(1, 0.45))

# Panel B with separated legends for terrestrial vs marine bioregions
description_time_region_cont_plot_leg <-
  plot_grid(description_time_region_cont_plot + theme(legend.position = "none"),
            plot_grid(get_legend(terr_description_legend_plot),
                      get_legend(marine_description_legend_plot),
                      ncol = 2),
            nrow = 2, rel_heights = c(1, 0.16))

# Combine figures
species_accumulation_combined_figure <- plot_grid(
  discovery_accumulation_alt_plot_leg,
  description_time_region_cont_plot_leg,
  nrow = 2, labels = "AUTO", rel_heights = c(0.4, 1),
  label_size = 10, label_fontfamily = "Sans", label_fontface = "bold", label_colour = "black")

ggsave(species_accumulation_combined_figure,
       filename = file.path(fig.dir, "species_accumulation_combined_figure.pdf"),
       width = 180, height = 210, units = "mm", dpi = 600, bg = "white", cairo_pdf)
ggsave(species_accumulation_combined_figure,
       filename = file.path(fig.dir, "species_accumulation_combined_figure.png"),
       width = 180, height = 210, units = "mm", dpi = 600, bg = "white")

# ==================================================================================================
## FIG S3 - SPECIES ACCUMULATION BY REGION - TERRESTRIAL GROUPS
# ==================================================================================================

### Calculate curves for each region by taxon combo
terr_taxa_description_curves_region_df <- taxa_full_list_clean %>% 
  left_join(taxonInfo) %>%
  subset(Region %in% terr_regions) %>%
  subset(! FacetGroup %in% c("Marine")) %>%
  subset(Present == 1) %>%
  ddply(.variables = .(FacetGroup, Taxon, Region),
        .fun = function(x){
          generate_description_curves(x$Year,
                                      min_year = subset(taxonomic_summary_statistics, Taxon == x$Taxon[1])$minYear,
                                      max_year = subset(taxonomic_summary_statistics, Taxon == x$Taxon[1])$maxYear)
        }) %>%
  ddply(.variables = .(FacetGroup, Taxon, Region),
        .fun = function(x){
          x$S_scaled <- x$S / max(x$S)
          return(x)
        }) %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order)) %>%
  mutate(Region = factor(Region,
                         levels = terr_regions,
                         labels = terr_region_labels)) %>%
  mutate(FacetGroup = factor(FacetGroup,
                             levels = facet_groups[1:3],
                             labels = c("Vascular plants",
                                        "Terrestrial & freshwater\nvertebrates",
                                        "Terrestrial & freshwater\ninvertebrates"))) %>%
  mutate(Region = factor(Region, levels = bioregion_info$Label))

### Plot
terr_discovery_accumulation_region_plot <- ggplot() +
  geom_path(data = terr_taxa_description_curves_region_df,
            aes(y = S_scaled, x = t, group = Taxon, colour = Taxon, linetype = Taxon),
            linewidth = 0.5) +
  labs(x = "Year", y = "Proportion of species described") +
  facet_grid(Region ~ FacetGroup, scale = "free_y") +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(1750, 1800, 1850, 1900, 1950, 2000),
                     limits = c(1750, 2025)) +
  scale_colour_manual(values = taxon_colors) +
  scale_linetype_manual(values = c(rep(1:2, 2)[1:4],     # Plants
                                   rep(1:2, 3)[1:5],     # Verts
                                   rep(1:2, 6)[1:11])) + # Inverts
  theme(panel.background   = element_blank(),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.background    = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.text    = element_text(family = "Sans", size = 5, colour = "black"),
        axis.title   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.ticks   = element_line(colour = "black", linewidth = 0.25),
        legend.title = element_blank(),# legend.title = element_text(family = "Sans", size = 7, colour = "black", hjust = 0.5),
        legend.text  = element_text(family = "Sans", size = 5, colour = "black"),
        legend.title.position = "top",
        legend.key = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())

### Generate separate legend for each group
for(i in 1:length(facet_groups[1:3])){
  test_plot <- ggplot() +
    geom_path(aes(y = S_scaled,
                  x = t,
                  group    = Taxon,
                  linetype = Taxon,
                  colour   = Taxon), 
    data = subset(terr_taxa_description_curves_region_df,
                  FacetGroup == facet_groups[i])) +
    # guides(colour   = guide_legend(ncol = ifelse(grepl(facet_groups[i], pattern = "invertebrates"), 2, 1)),
    #        linetype = guide_legend(ncol = ifelse(grepl(facet_groups[i], pattern = "invertebrates"), 2, 1))) +
    guides(colour   = guide_legend(facet_groups[i]),
           linetype = guide_legend(facet_groups[i])) +
    scale_colour_manual(values = subset(taxonInfo, FacetGroup == facet_groups[i])$Colour, ) +
    scale_linetype_manual(values = c(rep(1:2, 2)[1:4],     # Plants
                                     rep(1:2, 3)[1:5],     # Verts
                                     rep(1:2, 6)[1:11])) + # Inverts
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(family = "Arial", size = 7, colour = "black", hjust = 0.5),
          legend.text  = element_text(family = "Arial", size = 6, colour = "black"),
          legend.key.spacing.y = unit(0.01, "mm"),
          legend.box.margin = margin(c(0, 0, 0, 0)),
          legend.margin = margin(c(0, 0, 0, 0)),
          legend.key.size = unit(3, "mm"),
          legend.justification = "top",
          legend.justification.left = 1,
          legend.position = "right")
  FacetGroupLeg <- get_legend(test_plot)
  assign(paste("FacetGroupLeg", facet_group_abbrev[i], sep = "_"), FacetGroupLeg)
}

### Put plots and legend together
terr_discovery_accumulation_region_plot_leg <-
  plot_grid(terr_discovery_accumulation_region_plot + theme(legend.position = "none"),
            plot_grid(NULL,
                      FacetGroupLeg_plants, 
                      FacetGroupLeg_terr_fw_vertebrates,
                      FacetGroupLeg_terr_fw_invertebrates,
                      ncol = 1, align = "h", rel_heights = c(5, 4, 5, 11)),
            ncol = 2, rel_widths = c(1, 0.2))
terr_discovery_accumulation_region_plot_leg

ggsave(terr_discovery_accumulation_region_plot_leg,
       filename = file.path(fig.dir, "terr_discovery_accumulation_region.png"),
       height = 235, width = 180, units = "mm", dpi = 600, bg = "white")

ggsave(terr_discovery_accumulation_region_plot_leg,
       filename = file.path(fig.dir, "terr_discovery_accumulation_region.pdf"),
       height = 235, width = 180, units = "mm", dpi = 600, bg = "white", cairo_pdf)

# ==================================================================================================
## FIG S4 - SPECIES ACCUMULATION BY REGION - MARINE GROUPS
# ==================================================================================================

### Marine groups
marine_taxa_description_curves_region_df <- taxa_full_list_clean %>% 
  subset(Group %in% c("Marine")) %>%
  left_join(taxonInfo) %>%
  subset(Region %in% marine_regions) %>%
  subset(Present == 1) %>%
  ddply(.variables = .(FacetGroup, Taxon, Region),
        .fun = function(x){
          generate_description_curves(x$Year,
                                      min_year = subset(taxonomic_summary_statistics, Taxon == x$Taxon[1])$minYear,
                                      max_year = subset(taxonomic_summary_statistics, Taxon == x$Taxon[1])$maxYear)
        }) %>%
  ddply(.variables = .(FacetGroup, Taxon, Region),
        .fun = function(x){
          x$S_scaled <- x$S / max(x$S)
          return(x)
        }) %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order)) %>%
  mutate(Region = factor(Region,
                         levels = marine_regions,
                         labels = marine_region_labels))

### Plot
marine_discovery_accumulation_region_plot <- ggplot() +
  geom_path(aes(y = S_scaled, x = t, group= Taxon, colour = Taxon, linetype = Taxon),
            data = marine_taxa_description_curves_region_df,
            linewidth = 0.5) +
  labs(x = "Year", y = "Proportion of species described") +
  # facet_grid(Region ~ ., scale = "free_y") +
  facet_wrap(Region ~ ., scale = "free_y", ncol = 2) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(1750, 1800, 1850, 1900, 1950, 2000),
                     limits = c(1750, 2025)) +
  scale_colour_manual(values = taxon_colors[21:24]) +
  scale_linetype_manual(values = c(rep(1:2, 2)[1:4])) + # Marine
  theme(panel.background   = element_blank(),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.background    = element_blank(),
        panel.grid.major.x = element_line(color = "grey95"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.text    = element_text(family = "Sans", size = 5, colour = "black"),
        axis.title   = element_text(family = "Sans", size = 7, colour = "black"),
        axis.ticks   = element_line(colour = "black", linewidth = 0.25),
        legend.title = element_blank(),# legend.title = element_text(family = "Sans", size = 7, colour = "black", hjust = 0.5),
        legend.text  = element_text(family = "Sans", size = 5, colour = "black"),
        legend.title.position = "top",
        legend.key = element_blank(),
        legend.position = "right",
        strip.background = element_blank())
marine_discovery_accumulation_region_plot

ggsave(marine_discovery_accumulation_region_plot,
       filename = file.path(fig.dir, "marine_discovery_accumulation_region.png"),
       height = 180, width = 150, units = "mm", dpi = 600, bg = "white")

ggsave(marine_discovery_accumulation_region_plot,
       filename = file.path(fig.dir, "marine_discovery_accumulation_region.pdf"),
       # height = 235, width = 100, units = "mm", dpi = 600, bg = "white", cairo_pdf)
       height = 180, width = 180, units = "mm", dpi = 600, bg = "white", cairo_pdf)


# ==================================================================================================
## FIG. S5 - GENERATE NO OF DESCRIPTIONS PER YEAR PLOTS - SEPARATE TAXONOMIC GROUPS
# ==================================================================================================

taxa_temporal_counts <- taxa_full_list_clean %>% 
  select(all_of(c("Taxon","Species","Year"))) %>%
  subset(!is.na(Year)) %>%
  unique() %>%
  # ddply(.data = .,
  #       .variables = .(Taxon, Year),
  #       .fun = summarise,
  #       n_descriptions = length(unique(Species))) %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order))

description_time_cont_plot <- ggplot(data = taxa_temporal_counts) +
  # geom_bar(aes(x = Year, fill = Taxon)) +
  geom_bar(aes(x = Year, fill = Taxon), stat = "count", width = 1) +
  facet_wrap(~Taxon, nrow = 12, scale = "free_y") +
  scale_fill_manual(values = taxon_colors) +
  scale_x_continuous(breaks = seq(1750, 2020, 50),
                     limits = c(1750, 2025),
                     expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.1))) +
  labs(y = "No of species descriptions per year") +
  theme(panel.grid         = element_blank(),
        panel.background   = element_blank(),
        panel.border       = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        panel.grid.major.x = element_line(color = "grey95", linewidth = 0.5),
        panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.5),
        axis.ticks         = element_line(colour = "black", linewidth = 0.25),
        axis.text.x        = element_text(family = "Sans", size = 5, colour = "black", angle = 45, hjust = 1),
        axis.text.y        = element_text(family = "Sans", size = 5, colour = "black"),
        axis.title         = element_text(family = "Sans", size = 8, colour = "black"),
        strip.text         = element_text(family = "Sans", size = 7, colour = "black",
                                          margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "mm")),
        strip.background   = element_blank(),
        legend.position    = "none")
description_time_cont_plot

ggsave(description_time_cont_plot,
       filename = file.path(fig.dir, "description_time_cont.png"),
       height = 200, width = 180, units = "mm", dpi = 600, bg = "white")

ggsave(description_time_cont_plot,
       filename = file.path(fig.dir, "description_time_cont.pdf"),
       height = 200, width = 180, units = "mm", dpi = 600, bg = "white", cairo_pdf)

