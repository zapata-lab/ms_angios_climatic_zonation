# This files runs all the data wrangling and computationally intensive data 
# filtering for the manuscript 
#
# Executing this code generates `final_data_table.rds` (uncomment line 1020 to load), which 
# contains the data used in all analyses. In addition, it generates preliminary plots 
# (and even extra plots not used in the manuscript)
#
# By:   The authors
# Date: Aug 2024
#

#### Load libraries ####
library(tidyverse)
library(patchwork)
#library(car)
library(raster)
library(terra)
library(maps)
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2)
#install.packages("countrycode")
library(countrycode)
#install.packages("CoordinateCleaner")
library(CoordinateCleaner)
library(TreeTools)
#install.packages("rnaturalearthdata")
library(rnaturalearthdata)
#install.packages("geosphere")
library(geosphere)
library(sf)
library(diverge)
library(ggtree)

#### Preliminaries ####
time_start = Sys.time()
# Define `not in` to facilitate filtering
'%notin%' = Negate('%in%')

# Specify GBIF columns
columns_gbif = cols(
  gbifID = col_double(),
  datasetKey = col_character(),
  occurrenceID = col_character(),
  kingdom = col_character(),
  phylum = col_character(),
  class = col_character(),
  order = col_character(),
  family = col_character(),
  genus = col_character(),
  species = col_character(),
  infraspecificEpithet = col_skip(),
  taxonRank = col_skip(),
  scientificName = col_character(),
  verbatimScientificName = col_skip(),
  verbatimScientificNameAuthorship = col_skip(),
  countryCode = col_character(),
  locality = col_character(),
  stateProvince = col_character(),
  occurrenceStatus = col_skip(),
  individualCount = col_double(),
  publishingOrgKey = col_skip(),
  decimalLatitude = col_double(),
  decimalLongitude = col_double(),
  coordinateUncertaintyInMeters = col_double(),
  coordinatePrecision = col_double(),
  elevation = col_double(),
  elevationAccuracy = col_double(),
  depth = col_skip(),
  depthAccuracy = col_skip(),
  eventDate = col_datetime(format = ""),
  day = col_double(),
  month = col_double(),
  year = col_double(),
  taxonKey = col_double(),
  speciesKey = col_double(),
  basisOfRecord = col_character(),
  institutionCode = col_skip(),
  collectionCode = col_skip(),
  catalogNumber = col_skip(),
  recordNumber = col_skip(),
  identifiedBy = col_skip(),
  dateIdentified = col_skip(),
  license = col_skip(),
  rightsHolder = col_skip(),
  recordedBy = col_skip(),
  typeStatus = col_skip(),
  establishmentMeans = col_skip(),
  lastInterpreted = col_skip(),
  mediaType = col_skip(),
  issue = col_skip()
)

# Define colors for Figures
# Color-blind friendly palette http://www.cookbook-r.com/Graphs/Colors_(ggplot2)
pal = c("#000000", # Black
        "#E69F00", # Orange
        "#56B4E9", # Sky Blue
        "#009E73", # Bluish green
        "#F0E442", # Yellow
        "#0072B2", # Blue
        "#D55E00", # Vermillion
        "#CC79A7", # Reddish Purple
        "#999999") # Gray


#### Read in data ####
#----------------------------#
# # load gbif occurrences 
#----------------------------#
indata_monocots = read_tsv("../data//0194948-220831081235567.zip", 
                             col_types = columns_gbif)

indata_dicots = read_tsv("../data/0195061-220831081235567.zip", 
                           col_types = columns_gbif)

indata_angiosperms = rbind(indata_monocots,
                           indata_dicots)
#----------------------------#
# # load BioClim variables 
#----------------------------#
#

# To estimate Bioclimatic range
bio5 = raster("../data/wc2.1_30s_bio/wc2.1_30s_bio_5.tif") # Max Temperature of Warmest Month
bio6 = raster("../data/wc2.1_30s_bio/wc2.1_30s_bio_6.tif") # Min Temperature of Coldest Month
bio13 = raster("../data/wc2.1_30s_bio/wc2.1_30s_bio_13.tif") # Precipitation of Wettest Month
bio14 = raster("../data/wc2.1_30s_bio/wc2.1_30s_bio_14.tif") # Precipitation of Driest Month

# To estimate Bioclimatic overlap
# load Min Temp
tmin_01 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_01.tif") # Min temp Jan
tmin_02 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_02.tif") # Min temp Feb
tmin_03 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_03.tif") # Min temp Mar
tmin_04 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_04.tif") # Min temp Apr
tmin_05 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_05.tif") # Min temp May
tmin_06 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_06.tif") # Min temp Jun
tmin_07 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_07.tif") # Min temp Jul
tmin_08 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_08.tif") # Min temp Aug
tmin_09 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_09.tif") # Min temp Sep
tmin_10 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_10.tif") # Min temp Oct
tmin_11 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_11.tif") # Min temp Nov
tmin_12 = raster("../data/wc2.1_30s_tmin/wc2.1_30s_tmin_12.tif") # Min temp Dec
 
# load Max Temp
tmax_01 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_01.tif") # Max temp Jan
tmax_02 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_02.tif") # Max temp Feb
tmax_03 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_03.tif") # Max temp Mar
tmax_04 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_04.tif") # Max temp Apr
tmax_05 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_05.tif") # Max temp May
tmax_06 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_06.tif") # Max temp Jun
tmax_07 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_07.tif") # Max temp Jul
tmax_08 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_08.tif") # Max temp Aug
tmax_09 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_09.tif") # Max temp Sep
tmax_10 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_10.tif") # Max temp Oct
tmax_11 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_11.tif") # Min temp Nov
tmax_12 = raster("../data/wc2.1_30s_tmax/wc2.1_30s_tmax_12.tif") # Min temp Dec

#### Explore data ####
# Check tibbles
indata_monocots
indata_dicots
indata_angiosperms

# Plot occurrence data to get an overview
wm = borders("world", 
             colour = "gray50", 
             fill = "gray50")

# Monocots
# ggplot() + 
#   coord_fixed() + 
#   wm +
#   geom_point(data = indata_monocots, 
#              aes(x = decimalLongitude,
#                  y = decimalLatitude),
#              color = pal[2],
#              alpha = 0.3,
#              size = 0.5) +
#   theme_bw()
# 
# # Dicots
# ggplot() + 
#   coord_fixed() + 
#   wm +
#   geom_point(data = indata_dicots, 
#              aes(x = decimalLongitude, 
#                  y = decimalLatitude),
#              color = pal[6],
#              alpha = 0.3,
#              size = 0.5) +
#   theme_bw()

# Angiosperms
ggplot() + 
  coord_fixed() + 
  wm +
  geom_point(data = indata_angiosperms, 
             aes(x = decimalLongitude, 
                 y = decimalLatitude),
             color = pal[4],
             alpha = 0.3,
             size = 0.5) +
  theme_bw()


# Plot bioclimatic data.

# check the raster loaded correctly by plotting
plot(bio5, # Change to other variables
    main = "Bio5", # Change to other variables
    xlab = "Longitude",
    ylab = "Latitude",
    cex.axis = 1.3,
    cex.lab = 1.4,
    cex.main = 1.5,
    col = rev(heat.colors(10))
     )

#Check the raster loaded correctly by plotting
plot(tmax_12, # Change to other variables
     main = "tmax_12", # Change to other variables
     xlab = "Longitude",
     ylab = "Latitude",
     cex.axis = 1.3,
     cex.lab = 1.4,
     cex.main = 1.5,
     col = rev(heat.colors(10))
     )
 

#### Filter Data ####

#### Filter occurrence data ####

# Filter gbif species to species with at least 10 records (to
# try to avoid biases in range estimation)

indata_angiosperms = indata_angiosperms %>% 
  group_by(species) %>% 
  filter(n() >= 10)

# Confirm above filter worked
indata_angiosperms %>% 
  count(species) %>%
  arrange(n)

# Check that all records have coordinates
indata_angiosperms = indata_angiosperms %>%
  ungroup() %>%
  drop_na(decimalLongitude, 
          decimalLatitude,
          species,
          countryCode)

# Check that all records are HUMAN_OBSERVATION and PRESERVED_SPECIMEN.
# If not, apply filter
indata_angiosperms %>% 
  group_by(basisOfRecord) %>% 
  summarise(n = n())

# Check for individuals counts. Absence of records (count = 0) or really high
# records might indicate data entry problems
table(indata_angiosperms$individualCount)

# Remove individual count 0 and more than 500 (likely a mistake in data entry)
indata_angiosperms = indata_angiosperms %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 500 | is.na(individualCount))

# Remove data with too much coordinate uncertainty. USE 1Km MAXIMUM uncertainty.
# The coordinates are in meters, so if a coordinate has uncertainty 
# of 100 Km (100000 meters), divide uncertainty by 1000 (km) and filter by less than
# equal to 1 to use maximum uncertainty of 1km. Use 0.5 to filter for 500 m, 0.01 for 10 m
indata_angiosperms = indata_angiosperms %>%
  filter(coordinateUncertaintyInMeters/1000 <= 1 | # change for higher/lower uncertainty
           is.na(coordinateUncertaintyInMeters))

# To apply multiple filters implemented in `clean_coordinates` library,
# we need to make some edits to column names

# convert countryCode to countrycode (no capital letters)
names(indata_angiosperms)[12] = c("countrycode")


# convert country code from ISO2c to ISO3c
indata_angiosperms$countrycode = countrycode(indata_angiosperms$countrycode,
                                             origin =  "iso2c", 
                                             destination = "iso3c",
                                             nomatch = "XKX")

# check name changes
colnames(indata_angiosperms)

# convert Latitude and Longitude to latitude and longitude (no capital letters)
#names(indata_angiosperms)[16:17] = c("decimallatitude", 
#                                     "decimallongitude")

# Apply filters. For filter explanations, see https://ropensci.github.io/CoordinateCleaner/
indata_angiosperms = indata_angiosperms %>%
  cc_val() %>% # coordinate validity
  cc_equ() %>% # equal lat/long
  cc_cap(buffer = 10000) %>% # country capital using a buffer of 10Km around center of city
  cc_cen(buffer = 2000) %>% # country centroid using a buffer of 2Km from centroid
  cc_gbif(buffer = 2000) %>% # GBIF HQ using 2Km 
  cc_inst(buffer = 1000) %>% # Biodiversity Institution using 1Km
  cc_sea(scale = 110) %>% # restricting to land organisms. See https://ropensci.github.io/CoordinateCleaner/reference/cc_sea.html 
  cc_zero(buffer = 0.5) %>% # 0.5 decimal degrees around 0/0 lat/long
  cc_outl(method ="quantile", 
          mltpl = 5) %>% # outliers using quantile method and 5 times the Interquantile Range. See https://ropensci.github.io/CoordinateCleaner/reference/cc_outl.html
  cc_dupl() %>% # duplicates per species based on lat/long
  #cc_coun(iso3 = "countrycode") %>% # Lat/Long outside of country. With the update to this function, it broke. 
  cd_ddmm(ds = "datasetKey", 
          diff = 1) # error in degree conversion. See https://ropensci.github.io/CoordinateCleaner/reference/cd_ddmm.html

# apply filter based on elevation outliers and minimum sample size per species
indata_angiosperms = indata_angiosperms %>% 
  group_by(species) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) %>%
  mutate(across("class", ~str_replace(., "Liliopsida", "Monocots"))) %>%
  mutate(across("class", ~str_replace(., "Magnoliopsida", "Dicots"))) %>%
  dplyr::select(gbifID,
                class,
                family,
                genus,
                species,
                countrycode,
                decimalLatitude,
                decimalLongitude,
                elevation) %>%
ungroup()

# Up to this point, the table has 3,767,582 records

indata_angiosperms %>%
  count(species) %>%
  arrange(n)



#### Filter sister species pairs ####

# extract all sister pairs from the global tree in V.PhyloMaker2  
all_sp_pairs_phylo = extract_sisters(GBOTB.extended.TPL, sis_age = T)

# Make the previous dataframe a tibble and fix names to match in our dataset
all_sp_pairs_phylo = as_tibble(all_sp_pairs_phylo)
all_sp_pairs_phylo = all_sp_pairs_phylo %>%
  mutate(across("sp1", ~str_replace(., "_", " "))) %>%
  mutate(across("sp2", ~str_replace(., "_", " ")))


# Filter sister species pairs restricted to species of monocots with gbif data
#sp_pairs_monocots = all_sp_pairs_phylo %>%
#  filter(sp1 %in% indata_monocots$species & 
#           sp2 %in% indata_monocots$species)  

# Filter sister species pairs restricted to species of dicots with gbif data
#sp_pairs_dicots = all_sp_pairs_phylo %>%
#  filter(sp1 %in% indata_dicots$species & 
#           sp2 %in% indata_dicots$species)  

# join monocots and dicots
#sp_pairs_angiosperms = bind_rows(sp_pairs_monocots,
#                                 sp_pairs_dicots)



# Filter sister species pairs restricted to species with gbif data >= 10
sp_pairs_angiosperms = all_sp_pairs_phylo %>%
  filter(sp1 %in% indata_angiosperms$species & 
           sp2 %in% indata_angiosperms$species)  

# Make the list of pairs above longer and join all data
sp_pairs_angiosperms_long = pivot_longer(
  sp_pairs_angiosperms,
  cols = starts_with("sp"),
  names_to = "pair",
  values_to = "species") %>% 
  left_join(indata_angiosperms, 
            by = "species")  


#-----------------------------------------------#
#         Final table after all filters.        #
#-----------------------------------------------#

# Use this table to extract bioclimatic data for all species (in a pair)
sp_pairs_angiosperms_long = sp_pairs_angiosperms_long %>% 
  group_by(species)
sp_pairs_angiosperms_long
# 2734 species 684913 specimens


#### Filter bioclimatic data ####

# Make dataframes from coordinates to be able to extract values from raster layers
coords_all_sp_pairs_angiosperms = data.frame(
  lon = sp_pairs_angiosperms_long$decimalLongitude,
  lat = sp_pairs_angiosperms_long$decimalLatitude)

coordinates(coords_all_sp_pairs_angiosperms) = c("lon",
                                                 "lat")

# Plot to make sure it works
map()
points(coords_all_sp_pairs_angiosperms, 
       pch = 16, 
       col = pal[8])

# Extract values
bio5_val_all_sp_pairs_angiosperms = extract(x = bio5,
                                            y = coords_all_sp_pairs_angiosperms)
bio6_val_all_sp_pairs_angiosperms = extract(x = bio6,
                                            y = coords_all_sp_pairs_angiosperms)
bio13_val_all_sp_pairs_angiosperms = extract(x = bio13,
                                             y = coords_all_sp_pairs_angiosperms)
bio14_val_all_sp_pairs_angiosperms = extract(x = bio14,
                                             y = coords_all_sp_pairs_angiosperms)
tmin_01_val_all_sp_pairs_angiosperms = extract(x = tmin_01,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_02_val_all_sp_pairs_angiosperms = extract(x = tmin_02,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_03_val_all_sp_pairs_angiosperms = extract(x = tmin_03,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_04_val_all_sp_pairs_angiosperms = extract(x = tmin_04,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_05_val_all_sp_pairs_angiosperms = extract(x = tmin_05,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_06_val_all_sp_pairs_angiosperms = extract(x = tmin_06,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_07_val_all_sp_pairs_angiosperms = extract(x = tmin_07,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_08_val_all_sp_pairs_angiosperms = extract(x = tmin_08,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_09_val_all_sp_pairs_angiosperms = extract(x = tmin_09,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_10_val_all_sp_pairs_angiosperms = extract(x = tmin_10,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_11_val_all_sp_pairs_angiosperms = extract(x = tmin_11,
                                               y = coords_all_sp_pairs_angiosperms)
tmin_12_val_all_sp_pairs_angiosperms = extract(x = tmin_12,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_01_val_all_sp_pairs_angiosperms = extract(x = tmax_01,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_02_val_all_sp_pairs_angiosperms = extract(x = tmax_02,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_03_val_all_sp_pairs_angiosperms = extract(x = tmax_03,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_04_val_all_sp_pairs_angiosperms = extract(x = tmax_04,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_05_val_all_sp_pairs_angiosperms = extract(x = tmax_05,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_06_val_all_sp_pairs_angiosperms = extract(x = tmax_06,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_07_val_all_sp_pairs_angiosperms = extract(x = tmax_07,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_08_val_all_sp_pairs_angiosperms = extract(x = tmax_08,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_09_val_all_sp_pairs_angiosperms = extract(x = tmax_09,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_10_val_all_sp_pairs_angiosperms = extract(x = tmax_10,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_11_val_all_sp_pairs_angiosperms = extract(x = tmax_11,
                                               y = coords_all_sp_pairs_angiosperms)
tmax_12_val_all_sp_pairs_angiosperms = extract(x = tmax_12,
                                               y = coords_all_sp_pairs_angiosperms)

# Join with master tibble including all other data
sp_pairs_angiosperms_long_alldata = cbind(sp_pairs_angiosperms_long,
                                  bio5_val_all_sp_pairs_angiosperms,
                                  bio6_val_all_sp_pairs_angiosperms,
                                  bio13_val_all_sp_pairs_angiosperms,
                                  bio14_val_all_sp_pairs_angiosperms,
                                  tmin_01_val_all_sp_pairs_angiosperms,
                                  tmin_02_val_all_sp_pairs_angiosperms,
                                  tmin_03_val_all_sp_pairs_angiosperms,
                                  tmin_04_val_all_sp_pairs_angiosperms,
                                  tmin_05_val_all_sp_pairs_angiosperms,
                                  tmin_06_val_all_sp_pairs_angiosperms,
                                  tmin_07_val_all_sp_pairs_angiosperms,
                                  tmin_08_val_all_sp_pairs_angiosperms,
                                  tmin_09_val_all_sp_pairs_angiosperms,
                                  tmin_10_val_all_sp_pairs_angiosperms,
                                  tmin_11_val_all_sp_pairs_angiosperms,
                                  tmin_12_val_all_sp_pairs_angiosperms,
                                  tmax_01_val_all_sp_pairs_angiosperms,
                                  tmax_02_val_all_sp_pairs_angiosperms,
                                  tmax_03_val_all_sp_pairs_angiosperms,
                                  tmax_04_val_all_sp_pairs_angiosperms,
                                  tmax_05_val_all_sp_pairs_angiosperms,
                                  tmax_06_val_all_sp_pairs_angiosperms,
                                  tmax_07_val_all_sp_pairs_angiosperms,
                                  tmax_08_val_all_sp_pairs_angiosperms,
                                  tmax_09_val_all_sp_pairs_angiosperms,
                                  tmax_10_val_all_sp_pairs_angiosperms,
                                  tmax_11_val_all_sp_pairs_angiosperms,
                                  tmax_12_val_all_sp_pairs_angiosperms)

# Rename new columns 
sp_pairs_angiosperms_long_alldata  = sp_pairs_angiosperms_long_alldata  %>%
  rename(bio5 = ...13,
         bio6 = ...14,
         bio13 = ...15,
         bio14 = ...16,
         tmin01 = ...17,
         tmin02 = ...18,
         tmin03 = ...19,
         tmin04 = ...20,
         tmin05 = ...21,
         tmin06 = ...22,
         tmin07 = ...23,
         tmin08 = ...24,
         tmin09 = ...25,
         tmin10 = ...26,
         tmin11 = ...27,
         tmin12 = ...28,
         tmax01 = ...29,
         tmax02 = ...30,
         tmax03 = ...31,
         tmax04 = ...32,
         tmax05 = ...33,
         tmax06 = ...34,
         tmax07 = ...35,
         tmax08 = ...36,
         tmax09 = ...37,
         tmax10 = ...38,
         tmax11 = ...39,
         tmax12 = ...40)

# Remove missing data: localities without bioclimatic data
sp_pairs_angiosperms_long_alldata = sp_pairs_angiosperms_long_alldata %>%
  drop_na() 

# Check that for all species there are at least 10 specimens
sp_pairs_angiosperms_long_alldata %>%
  count(species) %>%
  arrange(n)


#### Elevation and Bioclimatic Ranges ####

#temporal objects to separate geographic regions
tmp_ntemp = sp_pairs_angiosperms_long_alldata %>% # This tibble is grouped by species. 
  ungroup() %>%
  group_by(tree_node) %>% # Important so sister pairs are not widespread across regions
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  group_by(species, tree_node, pair_age, pair, class) %>% # just for sorting final table and double check pairs are unique to region
  summarise(n_collections = n(), # count number of collections per species
            mean_latitude = mean(decimalLatitude), # a crude metric of geographic distribution
            elev_range = max(elevation) - min(elevation), # Elev range (response variable).
            temp_range = max(bio5) - min(bio6), # Temp range (response variable)
            precip_range = max(bio13) - min(bio14)) %>% # Ppt range (response variable)
  filter(elev_range > 0) %>% # in case there are species with 0 range
  add_column(region = "N. Temperate") %>% # assign region
  arrange(tree_node)

tmp_trop = sp_pairs_angiosperms_long_alldata %>% # This tibble is grouped by species.
  ungroup() %>%
  group_by(tree_node) %>% # Important so sister pairs are not widespread across regions
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics
  group_by(species, tree_node, pair_age, pair, class) %>% # just for sorting final table and double check pairs are unique to region
  summarise(n_collections = n(), # count number of collections per species
            mean_latitude = mean(decimalLatitude), # a crude metric of geographic distribution
            elev_range = max(elevation) - min(elevation), # Elev range (response variable).
            temp_range = max(bio5) - min(bio6), # Temp range (response variable)
            precip_range = max(bio13) - min(bio14)) %>% # Ppt range (response variable)
  filter(elev_range > 0) %>% # in case there are species with 0 range
  add_column(region = "Tropical") %>% # assign region
  arrange(tree_node)


tmp_stemp = sp_pairs_angiosperms_long_alldata %>% # This tibble is grouped by species.
  ungroup() %>%
  group_by(tree_node) %>% # Important so sister pairs are not widespread across regions
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  group_by(species, tree_node, pair_age, pair, class) %>% # just for sorting final table and double check pairs are unique to region
  summarise(n_collections = n(), # count number of collections per species
            mean_latitude = mean(decimalLatitude), # a crude metric of geographic distribution
            elev_range = max(elevation) - min(elevation), # Elev range (response variable).
            temp_range = max(bio5) - min(bio6), # Temp range (response variable)
            precip_range = max(bio13) - min(bio14)) %>% # Ppt range (response variable)
  filter(elev_range > 0) %>% # in case there are species with 0 range
  add_column(region = "S. Temperate") %>% # assign region
  arrange(tree_node)

# join all tibbles into a master tibble
sp_pairs_angiosperms_elev_biocl_range = bind_rows(
  tmp_ntemp,
  tmp_trop,
  tmp_stemp) %>% 
  ungroup()

# Check table for potential outliers
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = elev_range)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Elevation range") +
  theme(axis.text.x = element_text(size = 8)) 

# clearly a really strange plant with almost 7000 meters of range. Most likely an error in data
# entry. Detect outlier

sp_pairs_angiosperms_elev_biocl_range %>% 
  filter(elev_range == max(elev_range))

# remove outlier and sister using tree_node (NOTE Nov 2023: we already ran analyses with the outlier and results
# do not change )

sp_pairs_angiosperms_elev_biocl_range = sp_pairs_angiosperms_elev_biocl_range %>%
  filter(tree_node  != 102085)

# see table to count totals
sp_pairs_angiosperms_elev_biocl_range

# This table has 1486 species, N. Temp 1012, Tropics 190, and S. Temp 284
# Because the table with all data pre-filter by latitude has 2734 species, this means
# that 2734-1486 = 1248 species have sisters across latitudes

# uncomment below to confirm all species are unique to each region. n should start at 1
# sp_pairs_angiosperms_elev_biocl_range %>%
#  group_by(species) %>%
#  summarise(n = n()) %>%
#  arrange(desc(n))


# Table with all geographic data only for species restricted to each geographic region
sp_pairs_angiosperms_long_alldata_NotWidespread = sp_pairs_angiosperms_long_alldata %>% 
  filter(species %in% sp_pairs_angiosperms_elev_biocl_range$species)


# Summarize overall means per region
sp_pairs_angiosperms_elev_biocl_range_means_global = sp_pairs_angiosperms_elev_biocl_range %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev_range = mean(elev_range),
            sd_elev_range = sd(elev_range),
            mean_temp_range = mean(temp_range),
            sd_temp_range = sd(temp_range),
            mean_precip_range = mean(precip_range),
            sd_precip_range = sd(precip_range))

sp_pairs_angiosperms_elev_biocl_range_means_global

# Create table of counts and means per region and group
sp_pairs_angiosperms_elev_biocl_range_means_perGroup_perRegion = sp_pairs_angiosperms_elev_biocl_range %>%
  group_by(region, class) %>%
  summarise(n = n(),
            mean_elev_range = mean(elev_range),
            sd_elev_range = sd(elev_range),
            mean_temp_range = mean(temp_range),
            sd_temp_range = sd(temp_range),
            mean_precip_range = mean(precip_range),
            sd_precip_range = sd(precip_range))

sp_pairs_angiosperms_elev_biocl_range_means_perGroup_perRegion





#### Elevation and Bioclimatic Overlap ####

# Prepare data

tmp_summary_elev_tem_perSp = sp_pairs_angiosperms_long_alldata_NotWidespread %>% 
  summarise(min_elevation = min(elevation),
            max_elevation = max(elevation),
            mean_tmin01 = mean(tmin01), 
            mean_tmin02 = mean(tmin02),
            mean_tmin03 = mean(tmin03),
            mean_tmin04 = mean(tmin04),
            mean_tmin05 = mean(tmin05),
            mean_tmin06 = mean(tmin06),
            mean_tmin07 = mean(tmin07),
            mean_tmin08 = mean(tmin08),
            mean_tmin09 = mean(tmin09),
            mean_tmin10 = mean(tmin10),
            mean_tmin11 = mean(tmin11),
            mean_tmin12 = mean(tmin12),
            mean_tmax01 = mean(tmax01),
            mean_tmax02 = mean(tmax02),
            mean_tmax03 = mean(tmax03),
            mean_tmax04 = mean(tmax04),
            mean_tmax05 = mean(tmax05),
            mean_tmax06 = mean(tmax06),
            mean_tmax07 = mean(tmax07),
            mean_tmax08 = mean(tmax08),
            mean_tmax09 = mean(tmax09),
            mean_tmax10 = mean(tmax10),
            mean_tmax11 = mean(tmax11),
            mean_tmax12 = mean(tmax12),
            range_01 = mean_tmax01 - mean_tmin01,
            range_02 = mean_tmax02 - mean_tmin02,
            range_03 = mean_tmax03 - mean_tmin03,
            range_04 = mean_tmax04 - mean_tmin04,
            range_05 = mean_tmax05 - mean_tmin05,
            range_06 = mean_tmax06 - mean_tmin06,
            range_07 = mean_tmax07 - mean_tmin07,
            range_08 = mean_tmax08 - mean_tmin08,
            range_09 = mean_tmax09 - mean_tmin09,
            range_10 = mean_tmax10 - mean_tmin10,
            range_11 = mean_tmax11 - mean_tmin11,
            range_12 = mean_tmax12 - mean_tmin12)

# Join previous temp table with table summarizing elev and biocl range, region,
# age, pair, etc

sp_pairs_angiosperms_elev_biocl_overlap_perSp = left_join(sp_pairs_angiosperms_elev_biocl_range, 
                                                          tmp_summary_elev_tem_perSp, 
                                                          by = "species")

# Make the above table wider to see all pairs per row: USE THIS TABLE FOR ANALYSES
# In this table, the unique ID per row is tree_node (so, it's not included in the
# key:value pairs)

sp_pairs_angiosperms_elev_biocl_overlap_wide = pivot_wider(
  sp_pairs_angiosperms_elev_biocl_overlap_perSp,
  names_from = pair,
  values_from = c(species,
                  #tree_node,
                  pair_age,
                  class,
                  n_collections,
                  mean_latitude,
                  elev_range,
                  temp_range,
                  precip_range,
                  region,
                  min_elevation, 
                  max_elevation, 
                  mean_tmin01, 
                  mean_tmin02,
                  mean_tmin03,
                  mean_tmin04,
                  mean_tmin05,
                  mean_tmin06,
                  mean_tmin07,
                  mean_tmin08,
                  mean_tmin09,
                  mean_tmin10,
                  mean_tmin11,
                  mean_tmin12,
                  mean_tmax01,
                  mean_tmax02,
                  mean_tmax03,
                  mean_tmax04,
                  mean_tmax05,
                  mean_tmax06,
                  mean_tmax07,
                  mean_tmax08,
                  mean_tmax09,
                  mean_tmax10,
                  mean_tmax11,
                  mean_tmax12,
                  range_01,
                  range_02,
                  range_03,
                  range_04,
                  range_05,
                  range_06,
                  range_07,
                  range_08,
                  range_09,
                  range_10,
                  range_11,
                  range_12), 
  names_glue = "{pair}_{.value}")


# Estimate elevation and temperature overlap
# See Cadena et al 2011 PRSb and Kozak and Wiens 2007 PRSb for eqs.

sp_pairs_angiosperms_elev_biocl_range_overlap = sp_pairs_angiosperms_elev_biocl_overlap_wide %>% 
  rowwise() %>%
  mutate(elevation_overlap = (min(sp1_max_elevation, 
                                  sp2_max_elevation) - 
                                max(sp1_min_elevation, 
                                    sp2_min_elevation)) / 
           min(sp1_elev_range,
               sp2_elev_range),
         temp_overlap_jan = 0.5* ((min(sp1_mean_tmax01, 
                                       sp2_mean_tmax01) - 
                                     max(sp1_mean_tmin01, 
                                         sp2_mean_tmin01)) / 
                                    (sp1_range_01) +
                                    (min(sp1_mean_tmax01, 
                                         sp2_mean_tmax01) - 
                                       max(sp1_mean_tmin01, 
                                           sp2_mean_tmin01)) / 
                                    (sp2_range_01)),
         temp_overlap_feb = 0.5* ((min(sp1_mean_tmax02, 
                                       sp2_mean_tmax02) - 
                                     max(sp1_mean_tmin02, 
                                         sp2_mean_tmin02)) / 
                                    (sp1_range_02) +
                                    (min(sp1_mean_tmax02, 
                                         sp2_mean_tmax02) - 
                                       max(sp1_mean_tmin02, 
                                           sp2_mean_tmin02)) / 
                                    (sp2_range_02)),
         temp_overlap_mar = 0.5 * ((min(sp1_mean_tmax03, 
                                        sp2_mean_tmax03) -
                                      max(sp1_mean_tmin03, 
                                          sp2_mean_tmin03)) / 
                                     (sp1_range_03) +
                                     (min(sp1_mean_tmax03, 
                                          sp2_mean_tmax03) -
                                        max(sp1_mean_tmin03, 
                                            sp2_mean_tmin03)) /
                                     (sp2_range_03)),
         temp_overlap_apr = 0.5 * ((min(sp1_mean_tmax04, 
                                        sp2_mean_tmax04) - 
                                      max(sp1_mean_tmin04, 
                                          sp2_mean_tmin04)) / 
                                     (sp1_range_04) +
                                     (min(sp1_mean_tmax04, 
                                          sp2_mean_tmax04) -
                                        max(sp1_mean_tmin04, 
                                            sp2_mean_tmin04)) / 
                                     (sp2_range_04)),
         temp_overlap_may = 0.5 * ((min(sp1_mean_tmax05, 
                                        sp2_mean_tmax05) -
                                      max(sp1_mean_tmin05, 
                                          sp2_mean_tmin05)) / 
                                     (sp1_range_05) +
                                     (min(sp1_mean_tmax05,
                                          sp2_mean_tmax05) -
                                        max(sp1_mean_tmin05, 
                                            sp2_mean_tmin05)) / 
                                     (sp2_range_05)),
         temp_overlap_jun = 0.5 * ((min(sp1_mean_tmax06, 
                                        sp2_mean_tmax06) - 
                                      max(sp1_mean_tmin06, 
                                          sp2_mean_tmin06)) / 
                                     (sp1_range_06) +
                                     (min(sp1_mean_tmax06, 
                                          sp2_mean_tmax06) - 
                                        max(sp1_mean_tmin06, 
                                            sp2_mean_tmin06)) / 
                                     (sp2_range_06)),
         temp_overlap_jul = 0.5 * ((min(sp1_mean_tmax07, 
                                        sp2_mean_tmax07) - 
                                      max(sp1_mean_tmin07, 
                                          sp2_mean_tmin07)) / 
                                     (sp1_range_07) +
                                     (min(sp1_mean_tmax07, 
                                          sp2_mean_tmax07) -
                                        max(sp1_mean_tmin07, 
                                            sp2_mean_tmin07)) / 
                                     (sp2_range_07)),
         temp_overlap_aug = 0.5 * ((min(sp1_mean_tmax08,
                                        sp2_mean_tmax08) - 
                                      max(sp1_mean_tmin08, 
                                          sp2_mean_tmin08)) / 
                                     (sp1_range_08) +
                                     (min(sp1_mean_tmax08, 
                                          sp2_mean_tmax08) - 
                                        max(sp1_mean_tmin08, 
                                            sp2_mean_tmin08)) / 
                                     (sp2_range_08)),
         temp_overlap_sep = 0.5 * ((min(sp1_mean_tmax09, 
                                        sp2_mean_tmax09) -
                                      max(sp1_mean_tmin09, 
                                          sp2_mean_tmin09)) / 
                                     (sp1_range_09) +
                                     (min(sp1_mean_tmax09, 
                                          sp2_mean_tmax09) - 
                                        max(sp1_mean_tmin09, 
                                            sp2_mean_tmin09)) / 
                                     (sp2_range_09)),
         temp_overlap_oct = 0.5 * ((min(sp1_mean_tmax10, 
                                        sp2_mean_tmax10) -
                                      max(sp1_mean_tmin10, 
                                          sp2_mean_tmin10)) / 
                                     (sp1_range_10) +
                                     (min(sp1_mean_tmax10, 
                                          sp2_mean_tmax10) - 
                                        max(sp1_mean_tmin10, 
                                            sp2_mean_tmin10)) / 
                                     (sp2_range_10)),
         temp_overlap_nov = 0.5 * ((min(sp1_mean_tmax11,
                                        sp2_mean_tmax11) -
                                      max(sp1_mean_tmin11, 
                                          sp2_mean_tmin11)) / 
                                     (sp1_range_11) +
                                     (min(sp1_mean_tmax11, 
                                          sp2_mean_tmax11) - 
                                        max(sp1_mean_tmin11, 
                                            sp2_mean_tmin11)) / 
                                     (sp2_range_11)),
         temp_overlap_dec = 0.5 * ((min(sp1_mean_tmax12,
                                        sp2_mean_tmax12) - 
                                      max(sp1_mean_tmin12, 
                                          sp2_mean_tmin12)) / 
                                     (sp1_range_12) +
                                     (min(sp1_mean_tmax12, 
                                          sp2_mean_tmax12) - 
                                        max(sp1_mean_tmin12, 
                                            sp2_mean_tmin12)) / 
                                     (sp2_range_12))) %>%
  mutate(elevation_overlap = replace(elevation_overlap, 
                                     elevation_overlap < 0, 
                                     0),
         temp_overlap_jan = replace(temp_overlap_jan, 
                                    temp_overlap_jan < 0, 
                                    0),
         temp_overlap_feb = replace(temp_overlap_feb, 
                                    temp_overlap_feb < 0, 
                                    0),
         temp_overlap_mar = replace(temp_overlap_mar, 
                                    temp_overlap_mar < 0, 
                                    0),
         temp_overlap_apr = replace(temp_overlap_apr,
                                    temp_overlap_apr < 0,
                                    0),
         temp_overlap_may = replace(temp_overlap_may,
                                    temp_overlap_may < 0, 
                                    0),
         temp_overlap_jun = replace(temp_overlap_jun, 
                                    temp_overlap_jun < 0, 
                                    0),
         temp_overlap_jul = replace(temp_overlap_jul, 
                                    temp_overlap_jul < 0,
                                    0),
         temp_overlap_aug = replace(temp_overlap_aug, 
                                    temp_overlap_aug < 0, 
                                    0),
         temp_overlap_sep = replace(temp_overlap_sep, 
                                    temp_overlap_sep < 0, 
                                    0),
         temp_overlap_oct = replace(temp_overlap_oct, 
                                    temp_overlap_oct < 0, 
                                    0),
         temp_overlap_nov = replace(temp_overlap_nov, 
                                    temp_overlap_nov < 0, 
                                    0),
         temp_overlap_dec = replace(temp_overlap_dec, 
                                    temp_overlap_dec < 0, 
                                    0)) %>%
  mutate(temperature_overlap = sum(temp_overlap_jan,
                                   temp_overlap_feb,
                                   temp_overlap_mar,
                                   temp_overlap_apr,
                                   temp_overlap_may,
                                   temp_overlap_jun,
                                   temp_overlap_jul,
                                   temp_overlap_aug,
                                   temp_overlap_sep,
                                   temp_overlap_oct,
                                   temp_overlap_nov,
                                   temp_overlap_dec))

# THIS TABLE HAS ALL THE DATA WE NEED. IT'S ORGANIZED BY PAIRS IN EACH ROW 
# IT ALSO INCLUDES THE RESULTS FROM THE OVERLAP ANALYSES.
sp_pairs_angiosperms_elev_biocl_range_overlap
# 743 species pairs


# Visualize species chosen for analysis in the phylogeny
sp_included = sp_pairs_angiosperms_elev_biocl_range  %>%
  mutate(across("species", ~str_replace(., " ", "_")))


p = ggtree(GBOTB.extended.TPL, 
           size = 0.1,
           color = "gray70",
           layout = "circular")# +
  #geom_tiplab(align = TRUE)

p %<+% sp_included + 
  geom_tippoint(aes(color = region), 
                size = 1.5) +
  scale_colour_manual(values=c("#0072B2",
                               "#009E73",
                               "#D55E00"), 
                      na.translate=FALSE)


#............................................................................
# ------------------------------FINAL DATA TABLE-----------------------------
# ...........................................................................

# create a table just with the columns we need downstream for analysis and save it
final_data_table = sp_pairs_angiosperms_elev_biocl_range_overlap %>%
  dplyr::select(tree_node,
                sp1_pair_age,
                sp1_species,
                sp2_species,
                sp1_n_collections,
                sp2_n_collections,
                sp1_mean_latitude,
                sp2_mean_latitude,
                sp1_elev_range,
                sp2_elev_range,
                sp1_temp_range,
                sp2_temp_range,
                sp1_precip_range,
                sp2_precip_range,
                elevation_overlap,
                temperature_overlap,
                sp1_class,
                sp1_region) %>%
  rename(pair_age = sp1_pair_age,
         class = sp1_class,
         region = sp1_region) %>%
  ungroup()


# save the table above as final_data_table and load it again to facilitate analyses.
# Uncomment below to save
#saveRDS(final_data_table, "data/final_data_table.rds")

# uncomment below if want to load the data directly
#final_data_table = readRDS("data/final_data_table.rds")


# ...........................................................................
#...........................................................................
# ----------------------END CODE TO GENERATE DATA TABLE----------------------
# ...........................................................................
#...........................................................................




#### Figures: Range ####

# Boxplot elevation range per geographic region
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = elev_range)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Elevation range") +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot elevation range per geographic region and group
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = elev_range,
           col = class)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Elevation range") +
  scale_color_manual(values = c(pal[6], pal[2])) +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot temperature range per geographic region
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = temp_range)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Temperature range") +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot temperature range per geographic region and group
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = temp_range,
           col = class)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Temperature range") +
  scale_color_manual(values = c(pal[6], pal[2])) +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot precipitation range per geographic region
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = precip_range)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Precipitation range") +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot precipitation range per geographic region and group
ggplot(sp_pairs_angiosperms_elev_biocl_range,
       aes(x = region,
           y = precip_range,
           col = class)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Precipitation range") +
  scale_color_manual(values = c(pal[6], pal[2])) +
  theme(axis.text.x = element_text(size = 8)) 


#### Figures: Overlap ####

# Boxplot elevation overlap per geographic region
ggplot(final_data_table,
       aes(x = region,
           y = elevation_overlap,
           size = pair_age)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Elevation overlap") +
  theme(axis.text.x = element_text(size = 8))


# Boxplot elevation overlap per geographic region and group
ggplot(final_data_table,
       aes(x = region,
           y = elevation_overlap,
           col = class)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Elevation overlap") +
  scale_color_manual(values = c(pal[6], pal[2])) +
  theme(axis.text.x = element_text(size = 8)) 

# Boxplot temperature overlap per geographic region
ggplot(final_data_table,
       aes(x = region,
           y = temperature_overlap,
           size = pair_age)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Temperature overlap") +
  theme(axis.text.x = element_text(size = 8))


# Boxplot temperature overlap per geographic region and group
ggplot(final_data_table,
       aes(x = region,
           y = temperature_overlap,
           col = class)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_classic() +
  xlab("Geographic region") +
  ylab("Temperature overlap") +
  scale_color_manual(values = c(pal[6], pal[2])) +
  theme(axis.text.x = element_text(size = 8)) 


