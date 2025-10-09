# R scripts for Lim et al '*The rich, unique, and threatened biodiversity of tropical Asia*'

All scripts required for the analyses, results and figures in Lim et al '*The rich, unique, and threatened biodiversity of tropical Asia*'.

The project is organised into four folders which follow the broad sections of the manuscript:

1) Origins and assembly

2) Diversity and endemism

3) Biodiversity gaps

4) Threats

Within each folder is a ReadMe file outlining the scripts required. The order of running the scripts outlined in each ReadMe file is important - each script will process the raw data and produce intermediary data files necessary for generating the final figures.

Also, to run the 'Biodiversity gaps' and 'Threats' analyses, you will need to first generate the species intersections in the 'Diversity and endemism' folder.

IMPORTANT: most of the input data used in the analyses (such as IUCN range maps, GBIF data and checklist data sources), as well as some spatial data, such as GADM for coastlines, are required to run the scripts. As we do not own these data (and also many of them are very substantial - up to 250GB+) we can not provide them within this github project.

***To re-run the code the user will therefore have to obtain the data themselves and adjust the scripts as necessary to fit their personal folder structures and naming conventions***. These data sources are listed in the Supplementary Material, and also outlined in the ReadMe file of each of the four folders, plus at the top of the R scripts when necessary - all should be freely available from the respective data providers.
