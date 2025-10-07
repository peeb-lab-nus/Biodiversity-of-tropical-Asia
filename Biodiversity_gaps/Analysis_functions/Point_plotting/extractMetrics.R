#==================================================================================================#
#-------------------------------- Function for generating metrics ---------------------------------#
###
### Function for calculating sampling coverage metrics in 100 x 100 km grid cells
###
### Charlie Marsh
### charliem2003@github
###
#==================================================================================================#

extractMetrics <- function(pts) {
  taxonMap <- tibble()
  
  ################################################################################################
  ### Standardised sampling coverage - Total number of points per cell
  nPts <- pts %>%
    group_by(id_100km) %>%
    reframe(n_pts = n()) %>%
    mutate(log_pts = log(n_pts))
  
  ### Standardise by maximum no pts across all cells
  nPts <- nPts %>%
    mutate(value  = log_pts / max(log_pts)) %>%
    mutate(metric = "Standardised number of occurrences") %>%
    mutate(name   = taxon) %>%
    mutate(taxon  = taxonInfo$Taxon[taxonInfo$Name == taxon]) %>%
    mutate(group  = group_list[i]) %>%
    select(id_100km, value, metric, name, taxon, group)
  
  ### bind onto global dataframe
  taxonMap <- bind_rows(taxonMap, st_drop_geometry(nPts))
  
  ################################################################################################
  ### Taxonomic evenness - Proportion of singletons in each cell
  
  ### Only 1 record per year per species per cell
  sing <- pts %>%
    st_drop_geometry() %>%
    select(species, yearClean, id_100km) %>%
    distinct()
  
  ### Number of years with records for each species in each cell
  sing <- sing %>%
    group_by(id_100km, species) %>%
    reframe(n_pts = n())
  
  ### Calculate proportion of singletons
  sing <- sing %>%
    group_by(id_100km) %>%
    reframe(n_sing = sum(n_pts == 1),
            n_sp   = sum(n_pts > 0)) %>%
    mutate(prop_sing = (n_sing / n_sp)) # mutate(prop_sing = 1 - (n_sing / n_sp))
  
  ### Tidy up
  sing <- sing %>%
    mutate(value  = prop_sing) %>%
    mutate(metric = "Temporal replicates") %>%
    mutate(name   = taxon) %>%
    mutate(taxon  = taxonInfo$Taxon[taxonInfo$Name == taxon]) %>%
    mutate(group  = group_list[i]) %>%
    select(id_100km, value, metric, name, taxon, group)
  
  ### bind onto global dataframe
  taxonMap <- bind_rows(taxonMap, st_drop_geometry(sing))
  
  ################################################################################################
  ### Temporal coverage - Within each grid cell, find the most recent record for each species recorded
  
  ### Only 1 record per year per species per cell
  medYr <- pts %>%
    st_drop_geometry() %>%
    select(species, yearClean, id_100km) %>%
    distinct()
  
  ### For each cell the median record age
  medYr <- medYr %>%
    group_by(id_100km) %>%
    reframe(value = median(yearClean)) %>%
    mutate(metric = "Median record age") %>%
    mutate(name   = taxon) %>%
    mutate(taxon  = taxonInfo$Taxon[taxonInfo$Name == taxon]) %>%
    mutate(group  = group_list[i]) %>%
    select(id_100km, value, metric, name, taxon, group)
  
  ### bind onto global dataframe
  taxonMap <- bind_rows(taxonMap, st_drop_geometry(medYr))
}
