# PLOT ENDEMICITY ============

## PACKAGES ====================
rm(list = ls())

library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(sf)
library(cowplot)
library(patchwork)
library(ggnewscale)
library(scatterpie)

## DIRECTORIES ====================

# Locations of data, scripts and results - ADJUST FOR YOUR STRUCTURE
proj.dir <- "Diversity_and_endemism"

# You shouldn't need to adjust these folders
data.dir <- file.path(proj.dir, "Intersections", "Intersections")
res.dir  <- file.path(proj.dir, "Endemicity", "Results")
fig.dir  <- file.path(proj.dir, "Endemicity", "Figures")
if(!dir.exists(fig.dir)) { dir.create(fig.dir, recursive = TRUE) }

# ==================================================================================================
## DEFINE GRAPHICAL PARAMETERS
# ==================================================================================================

### Read in taxon info which has info on plot order and colours etc
taxon_info <- read.csv(file.path(proj.dir, "Global_totals.csv")) %>%
  arrange(desc(Group), desc(Taxon_Higher), Taxon) %>%
  mutate(Order = 1:nrow(.)) %>%
  mutate(Taxon = case_when(Taxon == "Flowering plants"       ~ "Flowering\nplants",
                           Taxon == "Freshwater fish"        ~ "Freshwater\nfish",
                           Taxon == "Freshwater crabs"       ~ "Freshwater\ncrabs",
                           Taxon == "Sharks and rays"        ~ "Sharks\nand rays",
                           Taxon == "Stick and leaf insects" ~ "Stick and\nleaf insects",
                           .default = Taxon)) %>%
  mutate(Taxon = factor(Taxon, levels = .$Taxon[.$Order]))

# Define factor orders and colours
taxon_order  <- taxon_info$Taxon
taxon_colors <- taxon_info$Colour

### Read in bioregion info which has info on plot order and colours etc
bioregion_info <- read.csv(file.path(proj.dir, "Bioregion_info.csv"))

### Define colours
group_colours <- taxon_info %>%
  mutate(Group = paste(Group, Taxon_Higher)) %>%
  group_by(Group) %>%
  reframe(col = rgb(apply(col2rgb(Colour), 1, mean)[1] / 255,
                    apply(col2rgb(Colour), 1, mean)[2] / 255,
                    apply(col2rgb(Colour), 1, mean)[3] / 255))
group_colours <- c("Plants" = group_colours$col[group_colours$Group == "Vascular plants "],
                   "Vert"   = group_colours$col[group_colours$Group == "Terrestrial and freshwater vertebrates"],
                   "Invert" = group_colours$col[group_colours$Group == "Terrestrial and freshwater invertebrates"],
                   "Marine" = group_colours$col[group_colours$Group == "Marine "])

# ==================================================================================================
## GENERATE DIVERSITY / ENDEMICITY METRICS
# ==================================================================================================

taxon_div_summary <- read.csv(file.path(res.dir, "taxon_div_summary.csv")) %>%
  mutate(Taxon = case_when(Taxon == "Flowering plants"       ~ "Flowering\nplants",
                           Taxon == "Freshwater fish"        ~ "Freshwater\nfish",
                           Taxon == "Freshwater crabs"       ~ "Freshwater\ncrabs",
                           Taxon == "Sharks and rays"        ~ "Sharks\nand rays",
                           Taxon == "Stick and leaf insects" ~ "Stick and\nleaf insects",
                           .default = Taxon)) %>%
  mutate(Taxon = factor(Taxon, levels = .$Taxon[taxon_info$Order]))

taxon_div_summary_final <- taxon_div_summary %>% 
  left_join(taxon_info) %>%
  mutate(n_native_nonend_species    = n_tropasia_species - n_tropasia_endemics) %>%
  mutate(prop_native_nonend_species = n_native_nonend_species / Global_total * 100) %>%
  mutate(prop_native_end_species    = n_tropasia_endemics     / Global_total * 100) %>%
  mutate(prop_global_div            = n_tropasia_species      / Global_total) %>%
  mutate(prop_tropasia_end          = n_tropasia_endemics     / n_tropasia_species)

### Split off terrestrial groups
terr_div_summary_final <- taxon_div_summary_final %>% 
  subset(Group != "Marine")

### Taxa with min and max proportion of global diversity
terr_div_summary_final$Taxon[which.min(terr_div_summary_final$prop_global_div)]
terr_div_summary_final$Taxon[which.max(terr_div_summary_final$prop_global_div)]

### Split off marine groups
mare_div_summary_final <- taxon_div_summary_final %>% 
  mutate(prop_global_div = n_tropasia_species / Global_total) %>%
  subset(Group == "Marine")

### Taxa with min and max proportion of global diversity
mare_div_summary_final$Taxon[which.min(mare_div_summary_final$prop_global_div)]
mare_div_summary_final$Taxon[which.max(mare_div_summary_final$prop_global_div)]

# ==================================================================================================
## PLOT GLOBAL ENDEMICITY
# ==================================================================================================

taxon_div_summary_final_melt <- taxon_div_summary_final %>%
  select("Taxon", "prop_native_nonend_species", "prop_native_end_species") %>%
  melt(id.vars = .(Taxon))

### Circular plot
taxon_order_padding <- c(as.character(taxon_info$Taxon[taxon_info$Group == "Vascular plants"]), NA,
                         as.character(taxon_info$Taxon[taxon_info$Taxon_Higher == "vertebrates"]), NA,
                         as.character(taxon_info$Taxon[taxon_info$Taxon_Higher == "invertebrates"]), NA,
                         as.character( taxon_info$Taxon[taxon_info$Group == "Marine"]), NA)
higher_grouping <- data.frame(group = factor(c("Plants", "Vert", "Invert", "Marine"),
                                             levels = c("Plants", "Vert", "Invert", "Marine")),
                              xmin = c(0, 5, 11, 23),
                              xmax = c(5, 11, 23, 28),
                              ymin = -70,
                              ymax = 80)

global_endemicity_plot <- ggplot() +
  theme(axis.text.theta = element_text(family = "Sans", size = 6, colour = "gray20",
                                       margin = margin(t = 0, r = 0, b = -5, 0, unit = "mm")),
        panel.background   = element_rect(fill = "white", color = "white"),
        plot.margin        = margin(t = 10, r = 10, b = 10, 10, unit = "mm"),
        panel.grid         = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title         = element_blank(),
        axis.ticks         = element_blank(),
        axis.text.y        = element_blank(),
        legend.position    = "none") +
  scale_x_discrete(limits = taxon_order_padding, breaks = taxon_order, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-70, 84), expand = c(0, 0)) +
  coord_radial(start = pi/180) +
  geom_hline(data = data.frame(y = c(0:4) * 20),
             aes(yintercept = y),
             color = "lightgrey") +
  scale_fill_manual(values = group_colours) +
  geom_rect(data = higher_grouping,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
            alpha = 0.3) +
  new_scale_fill() +
  scale_fill_manual(values = c(rep("grey95", length(taxon_order)), taxon_colors)) +
  geom_bar(data = taxon_div_summary_final_melt,
           aes(x = Taxon, y = value, fill = interaction(Taxon, variable)),
           stat = "identity", width = 0.9) +
  guides(theta = guide_axis_theta(angle = 90)) +
  annotate(x = 28, y = 24, label = "20%", geom = "text", color = "gray12", size = 2) +
  annotate(x = 28, y = 44, label = "40%", geom = "text", color = "gray12", size = 2) +
  annotate(x = 28, y = 64, label = "60%", geom = "text", color = "gray12", size = 2) +
  annotate(x = 28, y = 84, label = "80%", geom = "text", color = "gray12", size = 2) +
  annotate(x = 2.5, y = -25, label = "Vascular\nplants", geom = "text",
           family = "Sans", fontface = "bold", colour = "black", size = 2) +
  annotate(x = 8, y = -27, label = "Terrestrial &\nfreshwater\nvertebrates", geom = "text",
           family = "Sans", fontface = "bold", colour = "black", size = 2) +
  annotate(x = 17.5, y = -35, label = "Terrestrial &\nfreshwater\ninvertebrates", geom = "text",
           family = "Sans", fontface = "bold", colour = "black", size = 2) +
  annotate(x = 25.5, y = -25, label = "Marine", geom = "text",
           family = "Sans", fontface = "bold", colour = "black", size = 2)
global_endemicity_plot

ggsave(global_endemicity_plot,
       filename = file.path(fig.dir, "global_endemicity.pdf"),
       width = 120, height = 120, units = "mm", dpi = 600, bg = "white", Cairo::CairoPDF)
ggsave(global_endemicity_plot,
       filename = file.path(fig.dir, "global_endemicity.png"),
       width = 120, height = 120, units = "mm", dpi = 600, bg = "white")

# ==================================================================================================
## WITHIN-REGION ENDEMICITY SUMMARY 
# ==================================================================================================

# Read in maps
terr_bioregions <- readRDS(file.path(data.dir, "terr_bioregions.rds"))
mare_bioregions <- readRDS(file.path(data.dir, "mare_bioregions.rds"))
base_map        <- readRDS(file.path(data.dir, "gadm_base_map.rds"))

# Import data
region_end_invert <- read.csv(file.path(res.dir, "taxon_region_endemicity_long_df_terrestrial_and_freshwater_invertebrates_mean.csv"))
region_end_vert   <- read.csv(file.path(res.dir, "taxon_region_endemicity_long_df_terrestrial_and_freshwater_vertebrates_mean.csv"))
region_end_plant  <- read.csv(file.path(res.dir, "taxon_region_endemicity_long_df_vascular_plants_mean.csv"))
region_end_marine <- read.csv(file.path(res.dir, "taxon_region_endemicity_long_df_marine_mean.csv"))

region_endemism_all <- list(region_end_invert, region_end_vert, region_end_plant, region_end_marine) %>%
  do.call("rbind", .)

# Calculate area of bioregions
terr_area <- data.frame(X = terr_bioregions$X,
                        Y = terr_bioregions$Y,
                        Region = terr_bioregions$Bioregion,
                        AreaM2 = as.numeric(st_area(terr_bioregions)))
marine_area <- data.frame(X = mare_bioregions$X,
                          Y = mare_bioregions$Y,
                          Region = mare_bioregions$PROVINCE,
                          AreaM2 = as.numeric(st_area(mare_bioregions)))
combined_area <- rbind(terr_area, marine_area)

region_endemism_df <- left_join(region_endemism_all, combined_area, by = "Region") %>%
  mutate(LogArea = log(AreaM2 / 1000000)) %>%
  mutate(LogS = log(n_species)) %>%
  mutate(Endemic = n_endemics) %>%
  mutate(`Non endemic` = n_species - n_endemics) %>%
  mutate(FacetGroup = case_when(Group == "Marine" ~ "Marine",
                                Group == "Terrestrial and freshwater" & Taxon_Higher == "vertebrates"
                                ~ "Terrestrial &\nfreshwater vertebrates",
                                Group == "Terrestrial and freshwater" & Taxon_Higher == "invertebrates"
                                ~ "Terrestrial &\nfreshwater invertebrates",
                                Group == "Vascular plants" ~ "Vascular plants",
                                .default = NA))

# ==================================================================================================
## PLOT WITHIN-REGION ENDEMICITY MAPS (SUMMARY)
# ==================================================================================================

all_region_endemicity <- region_endemism_df %>%
  ddply(.data = .,
        .variables = .(FacetGroup, Region, X, Y),
        .fun = summarise, 
        Endemic = median(endemicity), 
        `Non endemic` = median(1 - endemicity))

facet_groups <- c("Vascular plants",
                  "Terrestrial &\nfreshwater vertebrates",
                  "Terrestrial &\nfreshwater invertebrates",
                  "Marine")
col_palettes <- c("Greens", "Purples", "YlOrBr", "Blues")

region_endemism_plot_list <- list()
for(i in 1:3){
  region_endemism_plot_list[[i]] <- ggplot() + 
    geom_sf(data = base_map,
            size = 0.1,
            color = NA,
            fill = "grey80") +
    geom_sf(aes(fill = Endemic),
            data = left_join(terr_bioregions, subset(all_region_endemicity, FacetGroup == facet_groups[i])), 
            color = NA) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_text(aes(x = 7714253, y = -1688714),
              label = facet_groups[i],
              hjust = 0, vjust = 0,
              color = "white", size.unit = "pt", size = 7, family = "Sans", fontface = "bold") +
    guides(fill = guide_colorbar(position = "inside",
                                 barwidth = 0.5,
                                 barheight = 2)) +
    scale_fill_distiller(name = "Endemicity",
                         palette = col_palettes[i],
                         limits = c(0, 0.84),
                         breaks = seq(0, 0.8, 0.2),
                         direction = 1) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(family = "Sans", colour = "white", size = 5, face = "bold"),
          legend.text = element_text(family = "Sans", colour = "white", size = 5, face = "bold"),
          legend.background = element_blank(),
          legend.frame = element_rect(colour = "white", linewidth = 0.1),
          legend.ticks = element_line(colour = "white", linewidth = 0.1),
          legend.ticks.length = unit(c(-1, 0), "mm"),
          legend.position.inside = c(0.999, 0.999),
          legend.justification.inside = c(1, 1),
          plot.margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm"),
          panel.background = element_rect(fill = "grey20")) +
    coord_sf(xlim = c(7614253, 17232260), ylim = c(-1788714, 3279685))
}

region_endemism_plot_list[[4]] <- ggplot() + 
  geom_sf(data = base_map,
          size = 0.1,
          fill = "grey20",
          color = NA) +
  geom_sf(aes(fill = Endemic),
          data = left_join(mare_bioregions,
                           subset(all_region_endemicity,
                                  FacetGroup == facet_groups[4]),
                           join_by(PROVINCE == Region)), 
          color = NA) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(x = 7714253, y = -1688714),
            label = facet_groups[4],
            hjust = 0, vjust = 0,
            color = "black", size.unit = "pt", size = 7, family = "Sans", fontface = "bold") +
  guides(fill = guide_colorbar(position = "inside",
                               barwidth = 0.5,
                               barheight = 2)) +
  scale_fill_distiller(name = "Endemicity",
                       palette = col_palettes[4],
                       limits = c(0, 0.2),
                       breaks = seq(0, 0.2, 0.05),
                       direction = 1) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(family = "Sans", colour = "black", size = 5, face = "bold"),
        legend.text = element_text(family = "Sans", colour = "black", size = 5, face = "bold"),
        legend.background = element_blank(),
        legend.frame = element_rect(colour = "black", linewidth = 0.1),
        legend.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.ticks.length = unit(c(-1, 0), "mm"),
        legend.position.inside = c(0.999, 0.999),
        legend.justification.inside = c(1, 1),
        plot.margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm"),
        panel.background = element_rect(fill = "grey80")) +
  coord_sf(xlim = c(7614253, 17232260), ylim = c(-1788714, 3279685))
region_endemism_plot_list[[4]]

region_endemism_combined_plot <- plot_grid(plotlist = region_endemism_plot_list,
                                           nrow = 2, align = "hv")
region_endemism_combined_plot

# ==================================================================================================
## FIGURE 3
# ==================================================================================================

fig3_combined <- global_endemicity_plot + region_endemism_combined_plot +
  plot_layout(nrow = 2,
              heights = unit(c(120, 70),  c("mm", "mm")),
              widths  = unit(120, "mm"),
              tag_level = 'new') +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.margin = margin(l = -10, t = 0, r = -10, b = 0, "mm"))) &
  theme(plot.tag.position = c(0.02, 1.00),
        plot.tag = element_text(family = "Sans", face = "bold", size = 10, hjust = 0, vjust = 0))
fig3_combined

ggsave(fig3_combined,
       filename = file.path(fig.dir, "fig3_combined.pdf"),
       width = 120, height = 200, units = "mm", dpi = 600, bg = "white", Cairo::CairoPDF)
ggsave(fig3_combined,
       filename = file.path(fig.dir, "fig3_combined.png"),
       width = 120, height = 200, units = "mm", dpi = 600, bg = "white")

# ==================================================================================================
## PLOT WITHIN-REGION ENDEMICITY MAPS
# ==================================================================================================

### Get dividing lines between subregions with land-sea border removed
base_lines <- st_union(st_cast(select(base_map, geom), "MULTILINESTRING"))
mare_lines <- st_union(st_cast(select(mare_bioregions, geom), "MULTILINESTRING"))
terr_lines <- st_union(st_cast(select(terr_bioregions, geom), "MULTILINESTRING"))

base_map_buffer <- st_buffer(base_lines, 10000)
terr_divides <- st_difference(terr_lines, base_map_buffer)
mare_divides <- st_difference(mare_lines, base_map_buffer)

### Subregion values
global_div <- taxon_info %>%
  mutate(Taxon = gsub("\n", " ", Taxon))
taxon_order <- global_div$Taxon

### Adjust position of Sumatra to avoid overlap
region_endemism_df <- region_endemism_df %>%
  mutate(X = case_when(Region == "Sumatra" ~ 11300000, .default = X)) %>%
  mutate(Y = case_when(Region == "Sumatra" ~  -250000, .default = Y))

endemicity_plot_list <- list()
for(i in 1:length(taxon_order)){
  if(subset(global_div, Taxon == taxon_order[i])$Group == "Marine"){
    # Mare
    base_plot <- ggplot() + 
      geom_sf(data = base_map, size = 0.1, fill = "grey20", color = NA) +
      geom_sf(data = mare_bioregions, color = NA) +
      geom_sf(data = mare_lines, col = "grey30", linewidth = 0.1) +
      guides(fill = "none") +
      scale_fill_distiller(palette = "PRGn") +
      theme(panel.background = element_rect(fill = "grey80"))
  } else {
    # Terr
    base_plot <- ggplot() + 
      geom_sf(data = base_map, size = 0.1, color = NA, fill = "grey80") +
      geom_sf(data = terr_bioregions, color = NA) +
      geom_sf(data = terr_lines, col = "grey30", linewidth = 0.1) +
      guides(fill = "none") +
      scale_fill_distiller(palette = "BrBG") +
      theme(panel.background = element_rect(fill = "white"))
  }
  endemicity_plot_list[[i]] <- base_plot +
    new_scale_fill() +
    geom_scatterpie(aes(x = X, y = Y, group = Region),
                    data = subset(region_endemism_df, Taxon == taxon_order[i]),
                    cols = c("Endemic", "Non endemic"),
                    color = NA, pie_scale = 2.5) +
    guides(fill = guide_legend(title = NULL)) +
    labs(title = taxon_order[i]) +
    scale_fill_manual(values = c("#003262","#FDB515")) +
    theme(panel.grid = element_blank(),
          # panel.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(family = "Sans", colour = "black", size = 8, hjust = 0.5, vjust = -1),
          plot.margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")
    ) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              col = "black", fill = NA) +
    coord_sf(xlim = c(7614253, 17232260),
             ylim = c(-1788714, 3400000),# 3279685),
             expand = FALSE)
}

### Export combined plot
combined_endemicity_plot_w_leg <-
  wrap_plots(endemicity_plot_list, ncol = 4, guides = "collect") &
  theme(legend.position = "bottom",
        legend.margin = margin(t = -1, b = 0, unit = "mm"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Sans", colour = "black", size = 8))
combined_endemicity_plot_w_leg

ggsave(combined_endemicity_plot_w_leg,
       filename = file.path(fig.dir, "endemicity_combined.pdf"),
       width = 180, height = 180, units = "mm", dpi = 600, bg = "white", Cairo::CairoPDF)
ggsave(combined_endemicity_plot_w_leg,
       filename = file.path(fig.dir, "endemicity_combined.png"),
       width = 180, height = 180, units = "mm", dpi = 600, bg = "white")
