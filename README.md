# R scripts for Lim et al '*Tropical Asia harbours a quarter of Earth’s species: knowledge gaps and priorities for protecting a global biodiversity hotspot*'

All scripts required for the analyses, results and figures in Lim et al '*Tropical Asia harbours a quarter of Earth’s species: knowledge gaps and priorities for protecting a global biodiversity hotspot*'.

## Instructions:

The project is organised into three folders which follow the broad sections of the manuscript:

1)  [Diversity and endemism](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/tree/main/Diversity_and_endemism)

2)  [Biodiversity gaps](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/tree/main/Biodiversity_gaps)

3)  [Threats](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/tree/main/Threats)

Within each folder is a ReadMe file outlining the scripts required. For detailed description of the analysis for each section, along with the data sources and locations required, refer to these ReadMe files for [Diversity and Endemism](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/tree/main/Diversity_and_endemism#readme), [Biodiversity Gaps](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/blob/main/Biodiversity_gaps/README.md), and [Threats](https://github.com/peeb-lab-nus/Biodiversity-of-tropical-Asia/blob/main/Threats/README.md).

The order of running the scripts outlined in each ReadMe file is important - each script will process the raw data and produce intermediary data files necessary for generating the final figures.

To run the 'Biodiversity gaps' and 'Threats' analyses, you will need to first generate the species intersections in the 'Diversity and endemism' folder.

IMPORTANT: most of the input data used in the analyses (such as IUCN range maps, GBIF data and checklist data sources), as well as some spatial data, such as GADM for coastlines, are required to run the scripts. As we do not own these data (and also many of them are very substantial - up to 250GB+) we can not provide them within this github project.

***To re-run the code the user will therefore have to obtain the data themselves from the sources outlined in Extended Data Table 1 and the Supplementary Material and adjust the scripts as necessary to fit their personal folder structures and naming conventions***. These data sources are listed in the Supplementary Material, and also outlined in the ReadMe file of each of the four folders, plus at the top of the R scripts when necessary - all should be freely available from the respective data providers. You will need the same versions of each data source, as outlined in the manuscript, to get the same results.

Some of the analysis steps, particularly those generating the intersections between species range maps and bioregions and the cleaning of GBIF data, are time-consuming and memory-intensive. Expect these steps to take a minimum of several weeks, and more likely months, depending on RAM - more RAM will allow for more parallelisation and will speed up the process.

## System requirements:

Code has been tested on Ubuntu 24.04, but some sections have also been tested on Windows and MacOS, and we expect all parts of the analysis should run on these OS systems.

No non-standard hardware is necessary, but RAM of 32+ GB is recommended for the larger spatial-based operations, and probably 128 GB is most likely necessary to clean the GBIF data. Significant storage space is required for the range maps, and particularly for the GBIF snapshot, which is \~260 GB.

## Installation guide:

A version of `R` 4.4+ is required (install from <https://cran.r-project.org/> for you OS).

The following R packages are necessary and can be installed through `install.packages()`: alphahull, arrow, Cairo, CoordinateCleaner, cowplot, dplyr, elevatr, ggnewscale, ggplot2, ggtext, httr2, iucnredlist, jsonlite, lubridate, parallel, patchwork, pbmcapply, plyr, purrr, readr, readxl, reshape2, scales, scatterpie, sf, shadowtext, sp, stringi, stringr, terra, tidyr, tidyterra, units.
