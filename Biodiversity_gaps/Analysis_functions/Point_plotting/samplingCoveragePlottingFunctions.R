#==================================================================================================#
#--------------------------------------- Plotting functions ---------------------------------------#
#==================================================================================================#

### Look up table for border colours depending on metric
borderCol <- function(metric) {
  if(metric == "C") {
    return("grey20")
  } else {
    return("grey80")
  }
}

### Function to plot one metric
occurrencePlot <- function(dd, facet) {
  ### Colour to fill in land and sea depends on taxa, and if it is citizen science-only data
  nrows <- nrow(unique(st_drop_geometry(dd[, facet])))
  
  if(facet == "taxon") {
    cols <- unique(st_drop_geometry(dd[, facet])) %>%
      left_join(st_drop_geometry(dd[, c("group", facet)]), keep = FALSE) %>%
      distinct() %>%
      mutate(colLand = case_when(group != "Marine" ~ "grey50",
                                 group == "Marine" ~ "grey90")) %>%
      mutate(colSea  = case_when(group != "Marine" ~ "white",
                                 group == "Marine" ~ "grey50"))
    colLand <- cols$colLand
    colSea  <- cols$colSea
  }
  if(facet == "group") {
    if(unique(dd$points) != "citizen") {
      colLand <- c(rep("grey50", nrows - 1), "grey90")
      colSea  <- c(rep("white",  nrows - 1), "grey50")
    } else {
      colLand <- c("grey50", "grey50", "grey50")
      colSea  <- c("white",  "white",  "white")
    }
  }
  
  p <- ggplot() +
    theme(panel.background      = element_blank(),
          panel.grid            = element_blank(),
          title                 = element_text(family = "Arial", colour = "black", face = "bold", size = 10, vjust = 0.5,
                                               margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "mm")
          ),
          strip.background      = element_blank(), #element_rect(fill = "white", colour = "black", linewidth = 0.2),
          strip.text            = element_text(family = "Arial", colour = "black", size = 8, vjust = 0,
                                               margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "mm")),
          panel.border          = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.margin           = margin(t = 0, r = 0, b = 0, l = -30, unit = "mm"),
          legend.position       = "bottom",
          legend.justification  = "center",
          legend.key.width      = unit(10, "mm"),
          legend.key.height     = unit(3, "mm"),
          legend.frame          = element_rect(colour = "black", linewidth = 0.2),
          legend.ticks          = element_line(colour = "black", linewidth = 0.2),
          legend.ticks.length   = unit(c(-1, 0), "mm"),
          legend.title.position = "bottom",
          legend.box.margin     = margin(t = -4, r = 0, b = 0, l = 0, unit = "mm"),
          legend.title          = element_text(family = "Arial", face = "plain", colour = "black", size = 6, hjust = 0.5,
                                               margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")),
          legend.text           = element_text(family = "Arial", colour = "black", size = 5),
          axis.text             = element_blank(), 
          axis.ticks            = element_blank()) +
    facet_wrap(. ~ get(facet),  ncol = 1, nrow = nrow(unique(st_drop_geometry(dd[, facet]))), strip.position = "left") +
    labs(title = unique(dd$metric)) +
    geom_sf(data = st_as_sfc(st_bbox(land)), colour = "black", fill = colSea) +
    geom_sf(data = land, colour = NA, linewidth = 0.05, fill = colLand) +
    scale_fill_viridis_c(option = "A") +
    geom_sf(data = dd, aes(fill = value), col = NA) +
    geom_sf(data = land, fill = NA, linewidth = 0.15, col = borderCol(unique(dd$metric))) +
    coord_sf(xlim = c(7614253, 17232260), ylim = c(-1788714, 3279685), expand = FALSE)
  return(p)
}

### Function to combine three metrics into final plot
panelPlot <- function(medsSplit, facet) {
  ### Loop through and plot each metric
  pAll <- purrr::map(medsSplit, ~ {
    occurrencePlot(dd = ., facet = facet)
  })
  
  ### Plot all together
  pAll <- patchwork::wrap_plots(pAll[[1]] +
                                  scale_fill_viridis_c(option = "A",
                                                       name = "Standardised number of occurrence records"),
                                pAll[[2]] +
                                  theme(strip.background = element_blank(),
                                        strip.text = element_blank()) + 
                                  scale_fill_viridis_c(option = "A", direction = -1,
                                                       name = "Prop. of temporal singletons"),
                                pAll[[3]] + 
                                  scale_fill_viridis_c(option = "A",
                                                       name = "Median record age") +
                                  theme(strip.background = element_blank(),
                                        strip.text = element_blank())) +
    plot_layout(guides = "keep")
  
  return(pAll)
}