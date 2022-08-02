# This kernel runs all the data wrangling and computationally intensive data filtering and analyses for the manuscript presented at https://github.com/zapata-lab/ms_angios_climatic_zonation
#
# Executing this code generates manuscript.RData, which contains analysis
# results. That file is then read by manuscript.rmd for rendering and
# presentation of the results.
#
# The code presented here is roughly in the order of the analyses presented in the manuscript, though there are exceptions.


#### Load libraries ####
library(tidyverse)
library(patchwork)
library(car)
library(raster)
library(maps)
library(V.PhyloMaker)
library(TreeTools)
library(rstatix)
library(diverge)
library(emmeans)

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
    individualCount = col_skip(),
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
    day = col_skip(),
    month = col_skip(),
    year = col_skip(),
    taxonKey = col_double(),
    speciesKey = col_double(),
    basisOfRecord = col_skip(),
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

# Load data with different Coordinate Uncertainties (CU) ####

# CU 10
all_monocots_CU10 = read_tsv("GBIF_DATA/All_Angios/ALL_Mono_10.zip", 
                           col_types = columns_gbif)
all_dicots_CU10 = read_tsv("GBIF_DATA/All_Angios/ALL_Dicot_10.zip", 
                         col_types = columns_gbif)

# CU 100
all_monocots_CU100 = read_tsv("GBIF_DATA/All_Angios/ALL_Mono_100.zip", 
                            col_types = columns_gbif)
all_dicots_CU100 = read_tsv("GBIF_DATA/All_Angios/ALL_Dicot_100.zip", 
                          col_types = columns_gbif)

# CU 500
all_monocots_CU500 = read_tsv("GBIF_DATA/All_Angios/ALL_Mono_500.zip", 
                            col_types = columns_gbif)
all_dicots_CU500 = read_tsv("GBIF_DATA/All_Angios/ALL_Dicot_500.zip", 
                          col_types = columns_gbif)

# CU 1000
all_monocots_CU1000 = read_tsv("GBIF_DATA/All_Angios/ALL_Mono_1000.zip", 
                             col_types = columns_gbif)
all_dicots_CU1000 = read_tsv("GBIF_DATA/All_Angios/ALL_Dicot_1000.zip", 
                           col_types = columns_gbif)


#### Data wrangling and filtering ####

# Monocots CU 10 ----
all_monocots_CU10_nomissID_unique_ptlsFiltered_tenPlus = all_monocots_CU10 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE)  %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 732 species, 68,220 records

# Dicots CU 10 ----
all_dicots_CU10_nomissID_unique_ptlsFiltered_tenPlus = all_dicots_CU10 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 2,397 species, 183,866 records

# Monocots CU 100 ----
all_monocots_CU100_nomissID_unique_ptlsFiltered_tenPlus = all_monocots_CU100 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 1,960 species, 276,033 records

# Dicots CU 100 ----
all_dicots_CU100_nomissID_unique_ptlsFiltered_tenPlus = all_dicots_CU100 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 7,783 species, 998,806 records

# Monocots CU 500 ----
all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus = all_monocots_CU500 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 2,266 species, 298,307 records

# Dicots CU 500 ----
all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus = all_dicots_CU500 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 9,248 species, 1,082,789 records

# Monocots CU 1000 ----
all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus = all_monocots_CU1000 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 2,936 species, 432,279 records

# Dicots CU 1000 ----
all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus = all_dicots_CU1000 %>%
  drop_na(species) %>% # drop NA
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 12,716 species, 1,545,480 records


# Analyses with CU 500 and CU 1000. Elevation Range ####

#------------------------------------------------#
# We are moving forward with CU 500 and CU 1000  #
#------------------------------------------------#

# Create temporal objects summarizing data per species and holding:
# - elevation range (our response variable of interest)
# - n_collections
# - taxonomic group
# - geographic region

# Monocots CU 500 ----

# N Temp ====
tmp_ntemp_monocots_CU500_perSp = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
  n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "N. Temperate") # assign region
  
# Resulting object has 646 species

# Tropics ====
tmp_trop_monocots_CU500_perSp = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 286 species

# S Temperate ====
tmp_stemp_monocots_CU500_perSp = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 675 species

# NOTE: the original object with all monocots had 2266 species. Separating the monocots by species in geographic regions results in: 646+286+675 = 1607. This means that 2266-1607 = 659 species were widespread across regions


# Dicots CU 500 ----

# N Temp ====
tmp_ntemp_dicots_CU500_perSp = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "N. Temperate") # assign region

# Resulting object has 2804 species

# Tropics ====
tmp_trop_dicots_CU500_perSp = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 1593 species

# S Temperate ====
tmp_stemp_dicots_CU500_perSp = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 2197 species

# NOTE: the original object with all dicots had 9248 species. Separating the dicots by species in geographic regions results in: 2804+1593+2197 = 6594. This mean that 9248-6594 = 2654 species were widespread across regions


# Join tibble of monocots

all_monocots_CU500_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU500_perSp, 
            tmp_trop_monocots_CU500_perSp,
            tmp_stemp_monocots_CU500_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_monocots_CU500_elevRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))

# Join tibble of dicots

all_dicots_CU500_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU500_perSp, 
            tmp_trop_dicots_CU500_perSp,
            tmp_stemp_dicots_CU500_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_dicots_CU500_elevRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))


# Join monocots and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU500_elevRange_coll_group_region = 
  bind_rows(all_monocots_CU500_elevRange_coll_group_region, 
            all_dicots_CU500_elevRange_coll_group_region)

# Resulting object has 8201 species


# Monocots CU 1000 ----

# N Temp ====
tmp_ntemp_monocots_CU1000_perSp = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "N. Temperate") # assign region

# Resulting object has 1188 species

# Tropics ====
tmp_trop_monocots_CU1000_perSp = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 322 species

# S Temperate ====
tmp_stemp_monocots_CU1000_perSp = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 694 species

# NOTE: the original object with all monocots had 2936 species. Separating the monocots by species in geographic regions results in: 1188+322+694  = 2204 This mean that 2936-2204 = 732 species were widespread across regions


# Dicots CU 1000 ----

# N Temp ====
tmp_ntemp_dicots_CU1000_perSp = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "N. Temperate") # assign region

# Resulting object has 5695 species

# Tropics ====
tmp_trop_dicots_CU1000_perSp = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 1864 species

# S Temperate ====
tmp_stemp_dicots_CU1000_perSp = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 2247 species

# NOTE: the original object with all dicots had 12716 species. Separating the dicots by species in geographic regions results in: 2247+5695+1864  = 9806 This mean that 12716-9806 = 2910 species were widespread across regions


# Join tibble of monocots

all_monocots_CU1000_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU1000_perSp, 
            tmp_trop_monocots_CU1000_perSp,
            tmp_stemp_monocots_CU1000_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_monocots_CU500_elevRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))

# Join tibble of dicots

all_dicots_CU1000_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU1000_perSp, 
            tmp_trop_dicots_CU1000_perSp,
            tmp_stemp_dicots_CU1000_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_dicots_CU500_elevRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))


# Join monocots and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU1000_elevRange_coll_group_region = 
  bind_rows(all_monocots_CU1000_elevRange_coll_group_region, 
            all_dicots_CU1000_elevRange_coll_group_region)

# resulting object has 12010 species



# Analysis Tropics south of Mexico ####

# Mexico is overcollected in the tropics; is the pattern in tropics driven by Mexico? Generate new tibbles for tropics with samples south of Mexico. Use CU 500 and CU 100

# South of Mexico

# Monocots CU 500 ----

# Tropical South of Mexico
tmp_trop_monocots_CU500_perSp_NoMex = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 14.533) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 12 species

# Dicots CU 500 ----

# Tropical south of Mexico
tmp_trop_dicots_CU500_perSp_NoMex = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 14.533) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 66 species

# Join tibble of monocots
all_monocots_CU500_NoMex_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU500_perSp, 
            tmp_trop_monocots_CU500_perSp_NoMex,
            tmp_stemp_monocots_CU500_perSp)

# Join tibble of dicots
all_dicots_CU500_NoMex_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU500_perSp, 
            tmp_trop_dicots_CU500_perSp_NoMex,
            tmp_stemp_dicots_CU500_perSp)


# Join monocts and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU500_NoMex_elevRange_coll_group_region = 
  bind_rows(all_monocots_CU500_NoMex_elevRange_coll_group_region, 
            all_dicots_CU500_NoMex_elevRange_coll_group_region)

# Resulting object has 6400 species

# summarise overall means per region
all_angiosperms_CU500_NoMex_elevRange_coll_group_region_counts_means_global = all_angiosperms_CU500_NoMex_elevRange_coll_group_region %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range))

# Create table of counts and means per region and group
all_angiosperms_CU500_NoMex_elevRange_coll_group_region_counts_means_perGroup_perRegion = all_angiosperms_CU500_NoMex_elevRange_coll_group_region %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range))


# Monocots CU 1000 ----

# Tropical South of Mexico
tmp_trop_monocots_CU1000_perSp_NoMex = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 14.533) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 31 species

# Dicots CU 1000 ----

# Tropical south of Mexico
tmp_trop_dicots_CU1000_perSp_NoMex = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 14.533) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 258 species

# Join tibble of monocots
all_monocots_CU1000_NoMex_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU1000_perSp, 
            tmp_trop_monocots_CU1000_perSp_NoMex,
            tmp_stemp_monocots_CU1000_perSp)

# Join tibble of dicots
all_dicots_CU1000_NoMex_elevRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU1000_perSp, 
            tmp_trop_dicots_CU1000_perSp_NoMex,
            tmp_stemp_dicots_CU1000_perSp)


# Join monocts and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU1000_NoMex_elevRange_coll_group_region = 
  bind_rows(all_monocots_CU1000_NoMex_elevRange_coll_group_region, 
            all_dicots_CU1000_NoMex_elevRange_coll_group_region)

# Resulting object has 10113 species

# summarise overall means per region
all_angiosperms_CU1000_NoMex_elevRange_coll_group_region_counts_means_global = all_angiosperms_CU1000_NoMex_elevRange_coll_group_region %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range))

# Create table of counts and means per region and group
all_angiosperms_CU1000_NoMex_elevRange_coll_group_region_counts_means_perGroup_perRegion = all_angiosperms_CU1000_NoMex_elevRange_coll_group_region %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range))

#-----------END ANALYSIS OF MEXICO--------------#


# Analysis considering strictly eudicots ====

non_eudicots = read_csv("orders_ANA_Magnoliids.csv")

# filter original dicots file to only include Eudicots

# Dicots CU 500
all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots = all_dicots_CU500 %>%
  drop_na(species) %>% # drop NA
  filter(order %notin% non_eudicots$order) %>% # filter non-eudicots
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 9,104 species, 1,054,174 records

# Dicots CU 1000
all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots = all_dicots_CU1000 %>%
  drop_na(species) %>% # drop NA
  filter(order %notin% non_eudicots$order) %>% # filter non eudicots
  group_by(species) %>% # group to select combination of:
  distinct(decimalLatitude, # unique lat
           decimalLongitude, # unique long
           eventDate, # unique dates. Important because records could be the same place but different collecting event
           .keep_all = TRUE) %>% 
  filter(elevation <= quantile(elevation, 0.975) & 
           elevation >= quantile(elevation, 0.025)) %>% # Filter elevation outliers (records in the 2.5% tails of the distribution per species)
  filter(n() >= 10) # filter poorly sampled species (species with less than 10 specimens)

# Resulting object: 12,545 species, 1,515,567 records



# Create temporal objects summarizing data per species and holding:
# - elevation range (our response variable of interest)
# - n_collections
# - taxonomic group
# - geographic region


# Eudicots CU 500 ----

# N Temp ====
tmp_ntemp_dicots_CU500_perSp_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "N. Temperate") # assign region

# Resulting object has 2798 species

# Tropics ====
tmp_trop_dicots_CU500_perSp_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 1538 species

# S Temperate ====
tmp_stemp_dicots_CU500_perSp_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 2145 species

# NOTE: the original object with all monocots had 9104 species. Separating the eudicots by species in geographic regions results in: 2798+1538+2145 = 6481. This mean that 9104-6481 = 2623 species were widespread across regions


# Join tibble of eudicots

all_dicots_CU500_elevRange_coll_group_region_onlyEudicots = 
  bind_rows(tmp_ntemp_dicots_CU500_perSp_onlyEudicots, 
            tmp_trop_dicots_CU500_perSp_onlyEudicots,
            tmp_stemp_dicots_CU500_perSp_onlyEudicots)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_dicots_CU500_elevRange_coll_group_region_onlyEudicots %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))


# Join monocots and eudicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots = 
  bind_rows(all_monocots_CU500_elevRange_coll_group_region, 
            all_dicots_CU500_elevRange_coll_group_region_onlyEudicots)

# Resulting object has 8088 species


# Eudicots CU 1000 ----

# N Temp ====
tmp_ntemp_dicots_CU1000_perSp_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) > 23.437) %>% # unique to N. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "N. Temperate") # assign region

# Resulting object has 5683 species

# Tropics ====
tmp_trop_dicots_CU1000_perSp_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437) %>% # unique to Tropics.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "Tropics") # assign region

# Resulting object has 1790 species

# S Temperate ====
tmp_stemp_dicots_CU1000_perSp_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(max(decimalLatitude) < -23.437) %>% # unique to S. Temp.
  summarise(elevation_range = max(elevation) - min(elevation), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  filter(elevation_range > 0) %>% # in case there are species with 0 range
  add_column(group = "Eudicots", # assign taxonomic group
             region = "S. Temperate") # assign region

# Resulting object has 2194 species

# NOTE: the original object with all eudicots had 12545 species. Separating the eudicots by species in geographic regions results in: 2194+1790+5683  = 9667 This mean that 12545-9667 = 2878 species were widespread across regions


# Join tibble of eudicots

all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots = 
  bind_rows(tmp_ntemp_dicots_CU1000_perSp_onlyEudicots, 
            tmp_trop_dicots_CU1000_perSp_onlyEudicots,
            tmp_stemp_dicots_CU1000_perSp_onlyEudicots)


# Join monocots and eudicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots = 
  bind_rows(all_monocots_CU1000_elevRange_coll_group_region, 
            all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots)

# resulting object has 11871 species

# We decided to run *all* analyses just with the eudicots, so the rest of the code for eudicots is presented below intermixed with the code for monocots and dicots



 # Bioclimatic analysis all species (Prediction 1) ####

#-------------------------------------------------------------------------------#
#
# BIOCLIM: Bioclimatic analysis
#
# Download the standard (19) WorldClim Bioclimatic variables for WorldClim version 2.1. These are the average for the years 1970-2000. Each download is a "zip" file containing 19 GeoTiff (.tif) files, one for each month of the variables. Taken from https://www.worldclim.org/data/worldclim21.html We downloaded the variables at 30 sec.
#
# Unzip the directory and use each raster to extract values
#
#-------------------------------------------------------------------------------#

# Get BioClim data and organize tables ====

#*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.#
#.*.*                                .*.*#
#.*.*         DO NOT RUN             .*.*#
#.*.*        THE CODE BELOW          .*.*#
#.*.*                                .*.*#
#*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.#


#The code below is *only for reference*. DO NOT RUN. To use the data, load the tables already with Bioclim variables at the very end; see below.

# Make tibbles with all records per species to extract Bioclim variables per record

# CU 500
# N Temp Monocots
# ntemp_monocots_CU500_bioclim = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#          species,
#          decimalLatitude,
#          decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "N. Temperate")
# 
# # Tropics Monocots
# trop_monocots_CU500_bioclim = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "Tropics")
# 
# # S Temp Monocots
# stemp_monocots_CU500_bioclim = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "S. Temperate")
# 
# # N. Temp Dicots
# ntemp_dicots_CU500_bioclim = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "N. Temperate")
# 
# ntemp_dicots_CU500_bioclim = ntemp_dicots_CU500_bioclim %>%
#   filter(species %in% tmp_ntemp_dicots_CU500_perSp$species)
# 
# # Tropics Dicots
# trop_dicots_CU500_bioclim = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "Tropics")
# 
# trop_dicots_CU500_bioclim = trop_dicots_CU500_bioclim %>%
#   filter(species %in% tmp_trop_dicots_CU500_perSp$species)
# 
# # S. Temperate Dicots
# stemp_dicots_CU500_bioclim = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "S. Temperate")
# 
# stemp_dicots_CU500_bioclim = stemp_dicots_CU500_bioclim %>%
#   filter(species %in% tmp_stemp_dicots_CU500_perSp$species)
# 
# 
# #----------#
# # EUDICOTS #
# #----------#
# 
# # N. Temp Eudicots
# ntemp_dicots_CU500_bioclim_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "N. Temperate")
# 
# # filter for only species present in region (remove widespread)
# ntemp_dicots_CU500_bioclim_onlyEudicots = ntemp_dicots_CU500_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_ntemp_dicots_CU500_perSp_onlyEudicots$species)
# 
# # Tropics Eudicots
# trop_dicots_CU500_bioclim_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "Tropics")
# 
# trop_dicots_CU500_bioclim_onlyEudicots = trop_dicots_CU500_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_trop_dicots_CU500_perSp_onlyEudicots$species)
# 
# # S. Temperate Eudicots
# stemp_dicots_CU500_bioclim_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "S. Temperate")
# 
# stemp_dicots_CU500_bioclim_onlyEudicots = stemp_dicots_CU500_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_stemp_dicots_CU500_perSp_onlyEudicots$species)
# 
# 
# # 
# # 
# # # CU 1000
# # # N Temp Monocots
# ntemp_monocots_CU1000_bioclim = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "N. Temperate")
# 
# # Tropics Monocots
# trop_monocots_CU1000_bioclim = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "Tropics")
# 
# # S Temp Monocots
# stemp_monocots_CU1000_bioclim = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Monocots",
#              region = "S. Temperate")
# 
# # N. Temp Dicots
# ntemp_dicots_CU1000_bioclim = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "N. Temperate")
# 
# ntemp_dicots_CU1000_bioclim = ntemp_dicots_CU1000_bioclim %>%
#   filter(species %in% tmp_ntemp_dicots_CU1000_perSp$species)
# 
# # Tropics Dicots
# trop_dicots_CU1000_bioclim = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "Tropics")
# 
# trop_dicots_CU1000_bioclim = trop_dicots_CU1000_bioclim %>%
#   filter(species %in% tmp_trop_dicots_CU1000_perSp$species)
# 
# # S. Temperate Dicots
# stemp_dicots_CU1000_bioclim = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Dicots",
#              region = "S. Temperate")
# 
# stemp_dicots_CU1000_bioclim = stemp_dicots_CU1000_bioclim %>%
#   filter(species %in% tmp_stemp_dicots_CU1000_perSp$species)
# 
# # 
# 
# #----------#
# # EUDICOTS #
# #----------#
# 
# 
# # N. Temp Dicots
# ntemp_dicots_CU1000_bioclim_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(min(decimalLatitude) > 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "N. Temperate")
# 
# ntemp_dicots_CU1000_bioclim_onlyEudicots = ntemp_dicots_CU1000_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_ntemp_dicots_CU1000_perSp_onlyEudicots$species)
# 
# # Tropics Eudicots
# trop_dicots_CU1000_bioclim_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(min(decimalLatitude) >= -23.437 &
#            max(decimalLatitude) <= 23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "Tropics")
# 
# trop_dicots_CU1000_bioclim_onlyEudicots = trop_dicots_CU1000_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_trop_dicots_CU1000_perSp_onlyEudicots$species)
# 
# # S. Temperate Eudicots
# stemp_dicots_CU1000_bioclim_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>%
#   filter(max(decimalLatitude) < -23.437) %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species) %>%
#   add_column(group = "Eudicots",
#              region = "S. Temperate")
# 
# stemp_dicots_CU1000_bioclim_onlyEudicots = stemp_dicots_CU1000_bioclim_onlyEudicots %>%
#   filter(species %in% tmp_stemp_dicots_CU1000_perSp_onlyEudicots$species)

#
#
#
#----------------------------#
# # load BioClim variables 
#----------------------------#
#
# bio5 = raster("wc2.1_30s_bio/wc2.1_30s_bio_5.tif") # Max Temperature of Warmest Month
# bio6 = raster("wc2.1_30s_bio/wc2.1_30s_bio_6.tif") # Min Temperature of Coldest Month
# bio13 = raster("wc2.1_30s_bio/wc2.1_30s_bio_13.tif") # Precipitation of Wettest Month
# bio14 = raster("wc2.1_30s_bio/wc2.1_30s_bio_14.tif") # Precipitation of Driest Month
# 
# # Check the raster loaded correctly by plotting
# # plot(bio5, # Change to other variables
# #     main = "Bio5", # Change to other variables
# #     xlab = "Longitude",
# #     ylab = "Latitude",
# #     cex.axis = 1.3,
# #     cex.lab = 1.4,
# #     cex.main = 1.5,
# #     col = rev(heat.colors(10))
# #     )
# 
# # N. Temp Monocots CU 500
# 
# # Make dataframes from coordinates to be able to extract values from raster layers
# coords_ntemp_monocots_CU500 = data.frame(
#   lon = ntemp_monocots_CU500_bioclim$decimalLongitude,
#   lat = ntemp_monocots_CU500_bioclim$decimalLatitude)
# coordinates(coords_ntemp_monocots_CU500) = c("lon",
#                                             "lat")
# # Plot to make sure it works
# #map()
# #points(coords_ntemp_monocots_CU500, pch=16)
# 
# # Extract values
# bio5_val_ntemp_monocots_CU500 = extract(x = bio5,
#                                        y = coords_ntemp_monocots_CU500)
# bio6_val_ntemp_monocots_CU500 = extract(x = bio6,
#                                         y = coords_ntemp_monocots_CU500)
# bio13_val_ntemp_monocots_CU500 = extract(x = bio13,
#                                         y = coords_ntemp_monocots_CU500)
# bio14_val_ntemp_monocots_CU500 = extract(x = bio14,
#                                         y = coords_ntemp_monocots_CU500)
# 
# # Join with tibble
# ntemp_monocots_CU500_bioclim = cbind(ntemp_monocots_CU500_bioclim,
#                                     bio5_val_ntemp_monocots_CU500,
#                                     bio6_val_ntemp_monocots_CU500,
#                                     bio13_val_ntemp_monocots_CU500,
#                                     bio14_val_ntemp_monocots_CU500)
# 
# ntemp_monocots_CU500_bioclim = ntemp_monocots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical Monocots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_monocots_CU500 = data.frame(
#   lon = trop_monocots_CU500_bioclim$decimalLongitude,
#   lat = trop_monocots_CU500_bioclim$decimalLatitude)
# coordinates(coords_trop_monocots_CU500) = c("lon",
#                                             "lat")
# # Extract values
# bio5_val_trop_monocots_CU500 = extract(x = bio5,
#                                        y = coords_trop_monocots_CU500)
# bio6_val_trop_monocots_CU500 = extract(x = bio6,
#                                         y = coords_trop_monocots_CU500)
# bio13_val_trop_monocots_CU500 = extract(x = bio13,
#                                          y = coords_trop_monocots_CU500)
# bio14_val_trop_monocots_CU500 = extract(x = bio14,
#                                          y = coords_trop_monocots_CU500)
# 
# # Join with tibble
# trop_monocots_CU500_bioclim = cbind(trop_monocots_CU500_bioclim,
#                                      bio5_val_trop_monocots_CU500,
#                                      bio6_val_trop_monocots_CU500,
#                                      bio13_val_trop_monocots_CU500,
#                                      bio14_val_trop_monocots_CU500)
# 
# trop_monocots_CU500_bioclim = trop_monocots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. Monocots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_monocots_CU500 = data.frame(
#   lon = stemp_monocots_CU500_bioclim$decimalLongitude,
#   lat = stemp_monocots_CU500_bioclim$decimalLatitude)
# coordinates(coords_stemp_monocots_CU500) = c("lon",
#                                             "lat")
# # Extract values
# bio5_val_stemp_monocots_CU500 = extract(x = bio5,
#                                        y = coords_stemp_monocots_CU500)
# bio6_val_stemp_monocots_CU500 = extract(x = bio6,
#                                        y = coords_stemp_monocots_CU500)
# bio13_val_stemp_monocots_CU500 = extract(x = bio13,
#                                         y = coords_stemp_monocots_CU500)
# bio14_val_stemp_monocots_CU500 = extract(x = bio14,
#                                         y = coords_stemp_monocots_CU500)
# 
# # Join with tibble
# stemp_monocots_CU500_bioclim = cbind(stemp_monocots_CU500_bioclim,
#                                     bio5_val_stemp_monocots_CU500,
#                                     bio6_val_stemp_monocots_CU500,
#                                     bio13_val_stemp_monocots_CU500,
#                                     bio14_val_stemp_monocots_CU500)
# 
# stemp_monocots_CU500_bioclim = stemp_monocots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # N. Temp Dicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_ntemp_dicots_CU500 = data.frame(
#   lon = ntemp_dicots_CU500_bioclim$decimalLongitude,
#   lat = ntemp_dicots_CU500_bioclim$decimalLatitude)
# coordinates(coords_ntemp_dicots_CU500) = c("lon",
#                                              "lat")
# 
# # Extract values
# bio5_val_ntemp_dicots_CU500 = extract(x = bio5,
#                                         y = coords_ntemp_dicots_CU500)
# bio6_val_ntemp_dicots_CU500 = extract(x = bio6,
#                                         y = coords_ntemp_dicots_CU500)
# bio13_val_ntemp_dicots_CU500 = extract(x = bio13,
#                                          y = coords_ntemp_dicots_CU500)
# bio14_val_ntemp_dicots_CU500 = extract(x = bio14,
#                                          y = coords_ntemp_dicots_CU500)
# 
# # Join with tibble
# ntemp_dicots_CU500_bioclim = cbind(ntemp_dicots_CU500_bioclim,
#                                      bio5_val_ntemp_dicots_CU500,
#                                      bio6_val_ntemp_dicots_CU500,
#                                      bio13_val_ntemp_dicots_CU500,
#                                      bio14_val_ntemp_dicots_CU500)
# 
# ntemp_dicots_CU500_bioclim = ntemp_dicots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical Dicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_dicots_CU500 = data.frame(
#   lon = trop_dicots_CU500_bioclim$decimalLongitude,
#   lat = trop_dicots_CU500_bioclim$decimalLatitude)
# coordinates(coords_trop_dicots_CU500) = c("lon",
#                                             "lat")
# # Extract values
# bio5_val_trop_dicots_CU500 = extract(x = bio5,
#                                        y = coords_trop_dicots_CU500)
# bio6_val_trop_dicots_CU500 = extract(x = bio6,
#                                        y = coords_trop_dicots_CU500)
# bio13_val_trop_dicots_CU500 = extract(x = bio13,
#                                         y = coords_trop_dicots_CU500)
# bio14_val_trop_dicots_CU500 = extract(x = bio14,
#                                         y = coords_trop_dicots_CU500)
# 
# # Join with tibble
# trop_dicots_CU500_bioclim = cbind(trop_dicots_CU500_bioclim,
#                                     bio5_val_trop_dicots_CU500,
#                                     bio6_val_trop_dicots_CU500,
#                                     bio13_val_trop_dicots_CU500,
#                                     bio14_val_trop_dicots_CU500)
# 
# trop_dicots_CU500_bioclim = trop_dicots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. Dicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_dicots_CU500 = data.frame(
#   lon = stemp_dicots_CU500_bioclim$decimalLongitude,
#   lat = stemp_dicots_CU500_bioclim$decimalLatitude)
# coordinates(coords_stemp_dicots_CU500) = c("lon",
#                                              "lat")
# # Extract values
# bio5_val_stemp_dicots_CU500 = extract(x = bio5,
#                                         y = coords_stemp_dicots_CU500)
# bio6_val_stemp_dicots_CU500 = extract(x = bio6,
#                                         y = coords_stemp_dicots_CU500)
# bio13_val_stemp_dicots_CU500 = extract(x = bio13,
#                                          y = coords_stemp_dicots_CU500)
# bio14_val_stemp_dicots_CU500 = extract(x = bio14,
#                                          y = coords_stemp_dicots_CU500)
# 
# # Join with tibble
# stemp_dicots_CU500_bioclim = cbind(stemp_dicots_CU500_bioclim,
#                                      bio5_val_stemp_dicots_CU500,
#                                      bio6_val_stemp_dicots_CU500,
#                                      bio13_val_stemp_dicots_CU500,
#                                      bio14_val_stemp_dicots_CU500)
# 
# stemp_dicots_CU500_bioclim = stemp_dicots_CU500_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 


#----------#
# EUDICOTS #
#----------#

# # N. Temp eudicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_ntemp_dicots_CU500_onlyEudicots = data.frame(
#   lon = ntemp_dicots_CU500_bioclim_onlyEudicots$decimalLongitude,
#   lat = ntemp_dicots_CU500_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_ntemp_dicots_CU500_onlyEudicots ) = c("lon",
#                                              "lat")
# 
# # Extract values
# bio5_val_ntemp_dicots_CU500_onlyEudicots  = extract(x = bio5,
#                                         y = coords_ntemp_dicots_CU500_onlyEudicots)
# bio6_val_ntemp_dicots_CU500_onlyEudicots  = extract(x = bio6,
#                                         y = coords_ntemp_dicots_CU500_onlyEudicots)
# bio13_val_ntemp_dicots_CU500_onlyEudicots  = extract(x = bio13,
#                                          y = coords_ntemp_dicots_CU500_onlyEudicots)
# bio14_val_ntemp_dicots_CU500_onlyEudicots  = extract(x = bio14,
#                                          y = coords_ntemp_dicots_CU500_onlyEudicots)
# 
# # Join with tibble
# ntemp_dicots_CU500_bioclim_onlyEudicots = cbind(ntemp_dicots_CU500_bioclim_onlyEudicots ,
#                                      bio5_val_ntemp_dicots_CU500_onlyEudicots,
#                                      bio6_val_ntemp_dicots_CU500_onlyEudicots,
#                                      bio13_val_ntemp_dicots_CU500_onlyEudicots,
#                                      bio14_val_ntemp_dicots_CU500_onlyEudicots)
# 
# ntemp_dicots_CU500_bioclim_onlyEudicots = ntemp_dicots_CU500_bioclim_onlyEudicots  %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical Eudicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_dicots_CU500_onlyEudicots = data.frame(
#   lon = trop_dicots_CU500_bioclim_onlyEudicots$decimalLongitude,
#   lat = trop_dicots_CU500_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_trop_dicots_CU500_onlyEudicots ) = c("lon",
#                                             "lat")
# # Extract values
# bio5_val_trop_dicots_CU500_onlyEudicots  = extract(x = bio5,
#                                        y = coords_trop_dicots_CU500_onlyEudicots )
# bio6_val_trop_dicots_CU500_onlyEudicots  = extract(x = bio6,
#                                        y = coords_trop_dicots_CU500_onlyEudicots )
# bio13_val_trop_dicots_CU500_onlyEudicots  = extract(x = bio13,
#                                         y = coords_trop_dicots_CU500_onlyEudicots )
# bio14_val_trop_dicots_CU500_onlyEudicots  = extract(x = bio14,
#                                         y = coords_trop_dicots_CU500_onlyEudicots )
# 
# # Join with tibble
# trop_dicots_CU500_bioclim_onlyEudicots = cbind(trop_dicots_CU500_bioclim_onlyEudicots ,
#                                     bio5_val_trop_dicots_CU500_onlyEudicots,
#                                     bio6_val_trop_dicots_CU500_onlyEudicots,
#                                     bio13_val_trop_dicots_CU500_onlyEudicots,
#                                     bio14_val_trop_dicots_CU500_onlyEudicots)
# 
# trop_dicots_CU500_bioclim_onlyEudicots = trop_dicots_CU500_bioclim_onlyEudicots  %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. eudicots CU 500
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_dicots_CU500_onlyEudicots = data.frame(
#   lon = stemp_dicots_CU500_bioclim_onlyEudicots$decimalLongitude,
#   lat = stemp_dicots_CU500_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_stemp_dicots_CU500_onlyEudicots ) = c("lon",
#                                              "lat")
# # Extract values
# bio5_val_stemp_dicots_CU500_onlyEudicots  = extract(x = bio5,
#                                         y = coords_stemp_dicots_CU500_onlyEudicots )
# bio6_val_stemp_dicots_CU500_onlyEudicots  = extract(x = bio6,
#                                         y = coords_stemp_dicots_CU500_onlyEudicots )
# bio13_val_stemp_dicots_CU500_onlyEudicots  = extract(x = bio13,
#                                          y = coords_stemp_dicots_CU500_onlyEudicots )
# bio14_val_stemp_dicots_CU500_onlyEudicots  = extract(x = bio14,
#                                          y = coords_stemp_dicots_CU500_onlyEudicots )
# 
# # Join with tibble
# stemp_dicots_CU500_bioclim_onlyEudicots = cbind(stemp_dicots_CU500_bioclim_onlyEudicots,
#                                      bio5_val_stemp_dicots_CU500_onlyEudicots,
#                                      bio6_val_stemp_dicots_CU500_onlyEudicots,
#                                      bio13_val_stemp_dicots_CU500_onlyEudicots,
#                                      bio14_val_stemp_dicots_CU500_onlyEudicots)
# 
# stemp_dicots_CU500_bioclim_onlyEudicots = stemp_dicots_CU500_bioclim_onlyEudicots %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)




# 
# # N. Temp Monocots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_ntemp_monocots_CU1000 = data.frame(
#   lon = ntemp_monocots_CU1000_bioclim$decimalLongitude,
#   lat = ntemp_monocots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_ntemp_monocots_CU1000) = c("lon",
#                                              "lat")
# 
# # Extract values
# bio5_val_ntemp_monocots_CU1000 = extract(x = bio5,
#                                         y = coords_ntemp_monocots_CU1000)
# bio6_val_ntemp_monocots_CU1000 = extract(x = bio6,
#                                         y = coords_ntemp_monocots_CU1000)
# bio13_val_ntemp_monocots_CU1000 = extract(x = bio13,
#                                          y = coords_ntemp_monocots_CU1000)
# bio14_val_ntemp_monocots_CU1000 = extract(x = bio14,
#                                          y = coords_ntemp_monocots_CU1000)
# 
# # Join with tibble
# ntemp_monocots_CU1000_bioclim = cbind(ntemp_monocots_CU1000_bioclim,
#                                      bio5_val_ntemp_monocots_CU1000,
#                                      bio6_val_ntemp_monocots_CU1000,
#                                      bio13_val_ntemp_monocots_CU1000,
#                                      bio14_val_ntemp_monocots_CU1000)
# 
# ntemp_monocots_CU1000_bioclim = ntemp_monocots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical Monocots CU 100
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_monocots_CU1000 = data.frame(
#   lon = trop_monocots_CU1000_bioclim$decimalLongitude,
#   lat = trop_monocots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_trop_monocots_CU1000) = c("lon",
#                                             "lat")
# # Extract values
# bio5_val_trop_monocots_CU1000 = extract(x = bio5,
#                                        y = coords_trop_monocots_CU1000)
# bio6_val_trop_monocots_CU1000 = extract(x = bio6,
#                                        y = coords_trop_monocots_CU1000)
# bio13_val_trop_monocots_CU1000 = extract(x = bio13,
#                                         y = coords_trop_monocots_CU1000)
# bio14_val_trop_monocots_CU1000 = extract(x = bio14,
#                                         y = coords_trop_monocots_CU1000)
# 
# # Join with tibble
# trop_monocots_CU1000_bioclim = cbind(trop_monocots_CU1000_bioclim,
#                                     bio5_val_trop_monocots_CU1000,
#                                     bio6_val_trop_monocots_CU1000,
#                                     bio13_val_trop_monocots_CU1000,
#                                     bio14_val_trop_monocots_CU1000)
# 
# trop_monocots_CU1000_bioclim = trop_monocots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. Monocots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_monocots_CU1000 = data.frame(
#   lon = stemp_monocots_CU1000_bioclim$decimalLongitude,
#   lat = stemp_monocots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_stemp_monocots_CU1000) = c("lon",
#                                              "lat")
# # Extract values
# bio5_val_stemp_monocots_CU1000 = extract(x = bio5,
#                                         y = coords_stemp_monocots_CU1000)
# bio6_val_stemp_monocots_CU1000 = extract(x = bio6,
#                                         y = coords_stemp_monocots_CU1000)
# bio13_val_stemp_monocots_CU1000 = extract(x = bio13,
#                                          y = coords_stemp_monocots_CU1000)
# bio14_val_stemp_monocots_CU1000 = extract(x = bio14,
#                                          y = coords_stemp_monocots_CU1000)
# 
# # Join with tibble
# stemp_monocots_CU1000_bioclim = cbind(stemp_monocots_CU1000_bioclim,
#                                      bio5_val_stemp_monocots_CU1000,
#                                      bio6_val_stemp_monocots_CU1000,
#                                      bio13_val_stemp_monocots_CU1000,
#                                      bio14_val_stemp_monocots_CU1000)
# 
# stemp_monocots_CU1000_bioclim = stemp_monocots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# 
# # N. Temp Dicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_ntemp_dicots_CU1000 = data.frame(
#   lon = ntemp_dicots_CU1000_bioclim$decimalLongitude,
#   lat = ntemp_dicots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_ntemp_dicots_CU1000) = c("lon",
#                                            "lat")
# 
# # Extract values
# bio5_val_ntemp_dicots_CU1000 = extract(x = bio5,
#                                       y = coords_ntemp_dicots_CU1000)
# bio6_val_ntemp_dicots_CU1000 = extract(x = bio6,
#                                       y = coords_ntemp_dicots_CU1000)
# bio13_val_ntemp_dicots_CU1000 = extract(x = bio13,
#                                        y = coords_ntemp_dicots_CU1000)
# bio14_val_ntemp_dicots_CU1000 = extract(x = bio14,
#                                        y = coords_ntemp_dicots_CU1000)
# 
# # Join with tibble
# ntemp_dicots_CU1000_bioclim = cbind(ntemp_dicots_CU1000_bioclim,
#                                    bio5_val_ntemp_dicots_CU1000,
#                                    bio6_val_ntemp_dicots_CU1000,
#                                    bio13_val_ntemp_dicots_CU1000,
#                                    bio14_val_ntemp_dicots_CU1000)
# 
# ntemp_dicots_CU1000_bioclim = ntemp_dicots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical Dicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_dicots_CU1000 = data.frame(
#   lon = trop_dicots_CU1000_bioclim$decimalLongitude,
#   lat = trop_dicots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_trop_dicots_CU1000) = c("lon",
#                                           "lat")
# # Extract values
# bio5_val_trop_dicots_CU1000 = extract(x = bio5,
#                                      y = coords_trop_dicots_CU1000)
# bio6_val_trop_dicots_CU1000 = extract(x = bio6,
#                                      y = coords_trop_dicots_CU1000)
# bio13_val_trop_dicots_CU1000 = extract(x = bio13,
#                                       y = coords_trop_dicots_CU1000)
# bio14_val_trop_dicots_CU1000 = extract(x = bio14,
#                                       y = coords_trop_dicots_CU1000)
# 
# # Join with tibble
# trop_dicots_CU1000_bioclim = cbind(trop_dicots_CU1000_bioclim,
#                                   bio5_val_trop_dicots_CU1000,
#                                   bio6_val_trop_dicots_CU1000,
#                                   bio13_val_trop_dicots_CU1000,
#                                   bio14_val_trop_dicots_CU1000)
# 
# trop_dicots_CU1000_bioclim = trop_dicots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. Dicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_dicots_CU1000= data.frame(
#   lon = stemp_dicots_CU1000_bioclim$decimalLongitude,
#   lat = stemp_dicots_CU1000_bioclim$decimalLatitude)
# coordinates(coords_stemp_dicots_CU1000) = c("lon",
#                                            "lat")
# # Extract values
# bio5_val_stemp_dicots_CU1000 = extract(x = bio5,
#                                       y = coords_stemp_dicots_CU1000)
# bio6_val_stemp_dicots_CU1000 = extract(x = bio6,
#                                       y = coords_stemp_dicots_CU1000)
# bio13_val_stemp_dicots_CU1000 = extract(x = bio13,
#                                        y = coords_stemp_dicots_CU1000)
# bio14_val_stemp_dicots_CU1000 = extract(x = bio14,
#                                        y = coords_stemp_dicots_CU1000)
# 
# # Join with tibble
# stemp_dicots_CU1000_bioclim = cbind(stemp_dicots_CU1000_bioclim,
#                                    bio5_val_stemp_dicots_CU1000,
#                                    bio6_val_stemp_dicots_CU1000,
#                                    bio13_val_stemp_dicots_CU1000,
#                                    bio14_val_stemp_dicots_CU1000)
# 
# stemp_dicots_CU1000_bioclim = stemp_dicots_CU1000_bioclim %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)



#----------#
# EUDICOTS #
#----------#

# # N. Temp Eudicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_ntemp_dicots_CU1000_onlyEudicots = data.frame(
#   lon = ntemp_dicots_CU1000_bioclim_onlyEudicots$decimalLongitude,
#   lat = ntemp_dicots_CU1000_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_ntemp_dicots_CU1000_onlyEudicots) = c("lon",
#                                            "lat")
# 
# # Extract values
# bio5_val_ntemp_dicots_CU1000_onlyEudicots = extract(x = bio5,
#                                       y = coords_ntemp_dicots_CU1000_onlyEudicots)
# bio6_val_ntemp_dicots_CU1000_onlyEudicots = extract(x = bio6,
#                                       y = coords_ntemp_dicots_CU1000_onlyEudicots)
# bio13_val_ntemp_dicots_CU1000_onlyEudicots = extract(x = bio13,
#                                        y = coords_ntemp_dicots_CU1000_onlyEudicots)
# bio14_val_ntemp_dicots_CU1000_onlyEudicots = extract(x = bio14,
#                                        y = coords_ntemp_dicots_CU1000_onlyEudicots)
# 
# # Join with tibble
# ntemp_dicots_CU1000_bioclim_onlyEudicots = cbind(ntemp_dicots_CU1000_bioclim_onlyEudicots,
#                                    bio5_val_ntemp_dicots_CU1000_onlyEudicots,
#                                    bio6_val_ntemp_dicots_CU1000_onlyEudicots,
#                                    bio13_val_ntemp_dicots_CU1000_onlyEudicots,
#                                    bio14_val_ntemp_dicots_CU1000_onlyEudicots)
# 
# ntemp_dicots_CU1000_bioclim_onlyEudicots = ntemp_dicots_CU1000_bioclim_onlyEudicots %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # Tropical eudicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_trop_dicots_CU1000_onlyEudicots = data.frame(
#   lon = trop_dicots_CU1000_bioclim_onlyEudicots$decimalLongitude,
#   lat = trop_dicots_CU1000_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_trop_dicots_CU1000_onlyEudicots) = c("lon",
#                                           "lat")
# # Extract values
# bio5_val_trop_dicots_CU1000_onlyEudicots = extract(x = bio5,
#                                      y = coords_trop_dicots_CU1000_onlyEudicots)
# bio6_val_trop_dicots_CU1000_onlyEudicots = extract(x = bio6,
#                                      y = coords_trop_dicots_CU1000_onlyEudicots)
# bio13_val_trop_dicots_CU1000_onlyEudicots = extract(x = bio13,
#                                       y = coords_trop_dicots_CU1000_onlyEudicots)
# bio14_val_trop_dicots_CU1000_onlyEudicots = extract(x = bio14,
#                                       y = coords_trop_dicots_CU1000_onlyEudicots)
# 
# # Join with tibble
# trop_dicots_CU1000_bioclim_onlyEudicots = cbind(trop_dicots_CU1000_bioclim_onlyEudicots,
#                                   bio5_val_trop_dicots_CU1000_onlyEudicots,
#                                   bio6_val_trop_dicots_CU1000_onlyEudicots,
#                                   bio13_val_trop_dicots_CU1000_onlyEudicots,
#                                   bio14_val_trop_dicots_CU1000_onlyEudicots)
# 
# trop_dicots_CU1000_bioclim_onlyEudicots = trop_dicots_CU1000_bioclim_onlyEudicots %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 
# 
# # S. Temp. eudicots CU 1000
# 
# # Make dataframes from coordinates to be able to use raster to extract values from raster layers
# coords_stemp_dicots_CU1000_onlyEudicots = data.frame(
#   lon = stemp_dicots_CU1000_bioclim_onlyEudicots$decimalLongitude,
#   lat = stemp_dicots_CU1000_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_stemp_dicots_CU1000_onlyEudicots) = c("lon",
#                                            "lat")
# # Extract values
# bio5_val_stemp_dicots_CU1000_onlyEudicots = extract(x = bio5,
#                                       y = coords_stemp_dicots_CU1000_onlyEudicots)
# bio6_val_stemp_dicots_CU1000_onlyEudicots = extract(x = bio6,
#                                       y = coords_stemp_dicots_CU1000_onlyEudicots)
# bio13_val_stemp_dicots_CU1000_onlyEudicots = extract(x = bio13,
#                                        y = coords_stemp_dicots_CU1000_onlyEudicots)
# bio14_val_stemp_dicots_CU1000_onlyEudicots = extract(x = bio14,
#                                        y = coords_stemp_dicots_CU1000_onlyEudicots)
# 
# # Join with tibble
# stemp_dicots_CU1000_bioclim_onlyEudicots = cbind(stemp_dicots_CU1000_bioclim_onlyEudicots,
#                                    bio5_val_stemp_dicots_CU1000_onlyEudicots,
#                                    bio6_val_stemp_dicots_CU1000_onlyEudicots,
#                                    bio13_val_stemp_dicots_CU1000_onlyEudicots,
#                                    bio14_val_stemp_dicots_CU1000_onlyEudicots)
# 
# stemp_dicots_CU1000_bioclim_onlyEudicots = stemp_dicots_CU1000_bioclim_onlyEudicots %>%
#   rename(bio5 = ...7,
#          bio6 = ...8,
#          bio13 = ...9,
#          bio14 = ...10)
# 




# Save tibbles as RDS
# saveRDS(stemp_dicots_CU1000_bioclim, "stemp_dicots_CU1000_bioclim.RDS")
# saveRDS(ntemp_dicots_CU1000_bioclim, "ntemp_dicots_CU1000_bioclim.RDS")
# saveRDS(trop_dicots_CU1000_bioclim, "trop_dicots_CU1000_bioclim.RDS")
# saveRDS(stemp_dicots_CU500_bioclim, "stemp_dicots_CU500_bioclim.RDS")
# saveRDS(ntemp_dicots_CU500_bioclim, "ntemp_dicots_CU500_bioclim.RDS")
# saveRDS(trop_dicots_CU500_bioclim, "trop_dicots_CU500_bioclim.RDS")
# saveRDS(stemp_monocots_CU500_bioclim, "stemp_monocots_CU500_bioclim.RDS")
# saveRDS(ntemp_monocots_CU500_bioclim, "ntemp_monocots_CU500_bioclim.RDS")
# saveRDS(trop_monocots_CU500_bioclim, "trop_monocots_CU500_bioclim.RDS")
# saveRDS(stemp_monocots_CU1000_bioclim, "stemp_monocots_CU1000_bioclim.RDS")
# saveRDS(ntemp_monocots_CU1000_bioclim, "ntemp_monocots_CU1000_bioclim.RDS")
# saveRDS(trop_monocots_CU1000_bioclim, "trop_monocots_CU1000_bioclim.RDS")

# saveRDS(stemp_dicots_CU1000_bioclim_onlyEudicots, "stemp_dicots_CU1000_bioclim_onlyEudicots.RDS")
# saveRDS(ntemp_dicots_CU1000_bioclim_onlyEudicots, "ntemp_dicots_CU1000_bioclim_onlyEudicots.RDS")
# saveRDS(trop_dicots_CU1000_bioclim_onlyEudicots, "trop_dicots_CU1000_bioclim_onlyEudicots.RDS")
# saveRDS(stemp_dicots_CU500_bioclim_onlyEudicots, "stemp_dicots_CU500_bioclim_onlyEudicots.RDS")
# saveRDS(ntemp_dicots_CU500_bioclim_onlyEudicots, "ntemp_dicots_CU500_bioclim_onlyEudicots.RDS")
# saveRDS(trop_dicots_CU500_bioclim_onlyEudicots, "trop_dicots_CU500_bioclim_onlyEudicots.RDS")




# Read in tibbles with Bioclim  data for all records. ----
# THESE ARE THE DATA THAT SHOULD BE USED TO ESTIMATE RANGES IN TEMPERATURE AND PRECIPITAION


#*****************************************************************#
#                                                                 #
#      YOU HAVE TO *RUN* THE FOLLOWING LINES TO LOAD THE DATA     #
#                                                                 #
#*****************************************************************#

ntemp_monocots_CU500_bioclim = readRDS("ntemp_monocots_CU500_bioclim.RDS")
trop_monocots_CU500_bioclim = readRDS("trop_monocots_CU500_bioclim.RDS")
stemp_monocots_CU500_bioclim = readRDS("stemp_monocots_CU500_bioclim.RDS")
ntemp_dicots_CU500_bioclim = readRDS("ntemp_dicots_CU500_bioclim.RDS")
trop_dicots_CU500_bioclim = readRDS( "trop_dicots_CU500_bioclim.RDS")
stemp_dicots_CU500_bioclim = readRDS("stemp_dicots_CU500_bioclim.RDS")
ntemp_monocots_CU1000_bioclim = readRDS("ntemp_monocots_CU1000_bioclim.RDS")
trop_monocots_CU1000_bioclim = readRDS("trop_monocots_CU1000_bioclim.RDS")
stemp_monocots_CU1000_bioclim = readRDS("stemp_monocots_CU1000_bioclim.RDS")
ntemp_dicots_CU1000_bioclim = readRDS("ntemp_dicots_CU1000_bioclim.RDS")
trop_dicots_CU1000_bioclim = readRDS("trop_dicots_CU1000_bioclim.RDS")
stemp_dicots_CU1000_bioclim = readRDS("stemp_dicots_CU1000_bioclim.RDS")

# Eudicots
stemp_dicots_CU1000_bioclim_onlyEudicots = readRDS("stemp_dicots_CU1000_bioclim_onlyEudicots.RDS")
ntemp_dicots_CU1000_bioclim_onlyEudicots = readRDS("ntemp_dicots_CU1000_bioclim_onlyEudicots.RDS")
trop_dicots_CU1000_bioclim_onlyEudicots = readRDS("trop_dicots_CU1000_bioclim_onlyEudicots.RDS")
stemp_dicots_CU500_bioclim_onlyEudicots = readRDS("stemp_dicots_CU500_bioclim_onlyEudicots.RDS")
ntemp_dicots_CU500_bioclim_onlyEudicots = readRDS("ntemp_dicots_CU500_bioclim_onlyEudicots.RDS")
trop_dicots_CU500_bioclim_onlyEudicots = readRDS("trop_dicots_CU500_bioclim_onlyEudicots.RDS")


# Estimate Thermal ranges ####

## Now estimate ranges in bioclimatic variables for each group, and then join tables as we did for elevation range

# Create temporal objects summarizing data per species and holding:
# - temp range (our response variable of interest)
# - precip range (our response variable of interest)
# - n_collections

# Monocots CU 500 ----

# N Temp ====
tmp_ntemp_monocots_CU500_bioclim_perSp = ntemp_monocots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 646 species

# Tropics ====
tmp_trop_monocots_CU500_bioclim_perSp = trop_monocots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 286 species

# S Temp ====
tmp_stemp_monocots_CU500_bioclim_perSp = stemp_monocots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 675 species


# Dicots CU 500 ----

# N Temp ====
tmp_ntemp_dicots_CU500_bioclim_perSp = ntemp_dicots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 2804 species

# Tropics ====
tmp_trop_dicots_CU500_bioclim_perSp = trop_dicots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 1593 species

# S Temp ====
tmp_stemp_dicots_CU500_bioclim_perSp = stemp_dicots_CU500_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 2197 species



## EUDICOTS 

# Eudicots CU 500 ====

# N Temp ====
tmp_ntemp_dicots_CU500_bioclim_perSp_onlyEudicots = ntemp_dicots_CU500_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 2798 species

# Tropics ====
tmp_trop_dicots_CU500_bioclim_perSp_onlyEudicots = trop_dicots_CU500_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 1538 species

# S Temp ====
tmp_stemp_dicots_CU500_bioclim_perSp_onlyEudicots = stemp_dicots_CU500_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 2145 species



# Join tibble of monocots

all_monocots_CU500_tempRange_precipRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU500_bioclim_perSp, 
            tmp_trop_monocots_CU500_bioclim_perSp,
            tmp_stemp_monocots_CU500_bioclim_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_monocots_CU500_tempRange_precipRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))

# Join tibble of dicots

all_dicots_CU500_tempRange_precipRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU500_bioclim_perSp, 
            tmp_trop_dicots_CU500_bioclim_perSp,
            tmp_stemp_dicots_CU500_bioclim_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_dicots_CU500_tempRange_precipRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))


# Join tibble of eudicots

all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots = 
  bind_rows(tmp_ntemp_dicots_CU500_bioclim_perSp_onlyEudicots, 
            tmp_trop_dicots_CU500_bioclim_perSp_onlyEudicots,
            tmp_stemp_dicots_CU500_bioclim_perSp_onlyEudicots)


# Join monocots and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU500_tempRange_precipRange_coll_group_region = 
  bind_rows(all_monocots_CU500_tempRange_precipRange_coll_group_region, 
            all_dicots_CU500_tempRange_precipRange_coll_group_region)

# Resulting object has 8201 species


# Join monocots and eudicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots = 
  bind_rows(all_monocots_CU500_tempRange_precipRange_coll_group_region, 
            all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)

# Resulting object has 8088 species



# Monocots CU 1000 ====

# N Temp ====
tmp_ntemp_monocots_CU1000_bioclim_perSp = ntemp_monocots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 1188 species

# Tropics ====
tmp_trop_monocots_CU1000_bioclim_perSp = trop_monocots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 322 species

# S Temp ====
tmp_stemp_monocots_CU1000_bioclim_perSp = stemp_monocots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Monocots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 694 species


# Dicots CU 1000 ====

# N Temp ====
tmp_ntemp_dicots_CU1000_bioclim_perSp = ntemp_dicots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 5695 species

# Tropics ====
tmp_trop_dicots_CU1000_bioclim_perSp = trop_dicots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 1864 species

# S Temp ====
tmp_stemp_dicots_CU1000_bioclim_perSp = stemp_dicots_CU1000_bioclim %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Dicots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 2247 species


## EUDICOTS

# Eudicots CU 1000 =====

# N Temp ====
tmp_ntemp_dicots_CU1000_bioclim_perSp_onlyEudicots = ntemp_dicots_CU1000_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "N. Temperate")

# Resulting object has 5683 species

# Tropics
tmp_trop_dicots_CU1000_bioclim_perSp_onlyEudicots = trop_dicots_CU1000_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "Tropics")

# Resulting object has 1790 species

# S Temp =====
tmp_stemp_dicots_CU1000_bioclim_perSp_onlyEudicots = stemp_dicots_CU1000_bioclim_onlyEudicots %>% # unique to N. Temp.
  drop_na() %>% 
  summarise(temp_range = max(bio5) - min(bio6), # estimate range, the response variable of interest.
            precip_range = max(bio13) - min(bio14), # estimate range, the response variable of interest.
            n_collections = n()) %>% # count number of collections per species
  add_column(group = "Eudicots", # assign taxonomic group
             region = "S. Temperate")

# Resulting object has 2194 species



# Join tibble of monocots

all_monocots_CU1000_tempRange_precipRange_coll_group_region = 
  bind_rows(tmp_ntemp_monocots_CU1000_bioclim_perSp, 
            tmp_trop_monocots_CU1000_bioclim_perSp,
            tmp_stemp_monocots_CU1000_bioclim_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_monocots_CCU100_tempRange_precipRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))

# Join tibble of dicots

all_dicots_CU1000_tempRange_precipRange_coll_group_region = 
  bind_rows(tmp_ntemp_dicots_CU1000_bioclim_perSp, 
            tmp_trop_dicots_CU1000_bioclim_perSp,
            tmp_stemp_dicots_CU1000_bioclim_perSp)

# Uncomment below to confirm records are unique to each region. n should be 1
#all_dicots_CU1000_tempRange_precipRange_coll_group_region %>% 
#  group_by(species) %>% 
#  summarise(n = n()) %>% 
#  arrange(desc(n))



# Join tibble of eudicots

all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots = 
  bind_rows(tmp_ntemp_dicots_CU1000_bioclim_perSp_onlyEudicots, 
            tmp_trop_dicots_CU1000_bioclim_perSp_onlyEudicots,
            tmp_stemp_dicots_CU1000_bioclim_perSp_onlyEudicots)


# Join monocots and dicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region = 
  bind_rows(all_monocots_CU1000_tempRange_precipRange_coll_group_region, 
            all_dicots_CU1000_tempRange_precipRange_coll_group_region)

# Resulting object has 12010 species


# Join monocts and eudicots to create final table of angiosperms
# ALL ANGIOSPERMS
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots = 
  bind_rows(all_monocots_CU1000_tempRange_precipRange_coll_group_region, 
            all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)

# Resulting object has 11871 species


# Create final tables with all data ####

#*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
# *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
#
#                      CREATE FINAL TABLES WITH *ALL* DATA
#
# *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
# *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 


# CU 500 ====
all_angiosperms_CU500_allData = all_angiosperms_CU500_elevRange_coll_group_region %>%
  left_join(all_angiosperms_CU500_tempRange_precipRange_coll_group_region, by = "species") %>%
  rename(n_collections = n_collections.x,
         group = group.x,
         region = region.x) %>%
  dplyr::select(species,
                group,
                region,
                n_collections,
                elevation_range,
                temp_range,
                precip_range)


# Summary tibble global
all_angiosperms_CU500_allData_counts_means_global = all_angiosperms_CU500_allData %>% group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))

# Summary tibble per group and per region
all_angiosperms_CU500_allData_counts_means_perGroup_perRegion = all_angiosperms_CU500_allData %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))


# CU 1000 ====
all_angiosperms_CU1000_allData = all_angiosperms_CU1000_elevRange_coll_group_region %>%
  left_join(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region, by = "species") %>%
  rename(n_collections = n_collections.x,
         group = group.x,
         region = region.x) %>%
  dplyr::select(species,
                group,
                region,
                n_collections,
                elevation_range,
                temp_range,
                precip_range)

# Summary tibble global
all_angiosperms_CU1000_allData_counts_means_global = all_angiosperms_CU1000_allData %>% group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))

# Summary tibble per group and per region
all_angiosperms_CU1000_allData_counts_means_perGroup_perRegion = all_angiosperms_CU1000_allData %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))



### ### ### ### ### ### ### 
### ONLY WITH EUDICOTS  ### 
### ### ### ### ### ### ### 

# CU 500 (Eudicots) ====
all_angiosperms_CU500_allData_onlyEudicots = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots %>%
  left_join(all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots, by = "species") %>%
  rename(n_collections = n_collections.x,
         group = group.x,
         region = region.x) %>%
  dplyr::select(species,
                group,
                region,
                n_collections,
                elevation_range,
                temp_range,
                precip_range)


# Summary tibble global
all_angiosperms_CU500_allData_counts_means_global_onlyEudicots = all_angiosperms_CU500_allData_onlyEudicots %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))

# Summary tibble per group and per region
all_angiosperms_CU500_allData_counts_means_perGroup_perRegion_onlyEudicots = all_angiosperms_CU500_allData_onlyEudicots %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))


# CU 1000 (Eudicots) ====
all_angiosperms_CU1000_allData_onlyEudicots = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots %>%
  left_join(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots, by = "species") %>%
  rename(n_collections = n_collections.x,
         group = group.x,
         region = region.x) %>%
  dplyr::select(species,
                group,
                region,
                n_collections,
                elevation_range,
                temp_range,
                precip_range)

# Summary tibble global
all_angiosperms_CU1000_allData_counts_means_global_onlyEudicots = all_angiosperms_CU1000_allData_onlyEudicots %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))

# Summary tibble per group and per region
all_angiosperms_CU1000_allData_counts_means_perGroup_perRegion_onlyEudicots = all_angiosperms_CU1000_allData_onlyEudicots %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range),
            sd_elev_range = sd(elevation_range),
            sd_temp_range = sd(temp_range),
            sd_precip_range = sd(precip_range))


# Extract genera for species-pairs analyses ####

#-------------------------------------------------------------------------#
#
# Extract unique genera to look for species pairs in the plant phylogeny
#
#-------------------------------------------------------------------------#

# CU 500 ####

## N Temp Monocots CU500
genera_ntemp_monocots_CU500 = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop Monocots CU500
genera_trop_monocots_CU500 = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp Monocots CU500
genera_stemp_monocots_CU500 = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## N Temp Dicots CU500
genera_ntemp_dicots_CU500 = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop Dicots CU500
genera_trop_dicots_CU500 = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp Dicots CU500
genera_stemp_dicots_CU500 = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)


## N Temp eudicots CU500 - only EUDICOTS
genera_ntemp_dicots_CU500_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop eudicots CU500 - only EUDICOTS
genera_trop_dicots_CU500_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp eudicots CU500 - only EUDICOTS
genera_stemp_dicots_CU500_onlyEudicots = all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)


# CU 1000 ####

## N Temp Monocots CU1000
genera_ntemp_monocots_CU1000 = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop Monocots CU1000
genera_trop_monocots_CU1000 = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp Monocots CU1000
genera_stemp_monocots_CU1000 = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## N Temp Dicots CU1000
genera_ntemp_dicots_CU1000 = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop Dicots CU1000
genera_trop_dicots_CU1000 = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp Dicots CU1000
genera_stemp_dicots_CU1000 = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)


## N Temp eudicots CU1000 - only EUDICOTS
genera_ntemp_dicots_CU1000_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) > 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## Trop eudicots CU1000 - only EUDICOTS
genera_trop_dicots_CU1000_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(min(decimalLatitude) >= -23.437 & 
           max(decimalLatitude) <= 23.437)  %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)

## S Temp Eudiicots CU1000 - only EUDICOTS
genera_stemp_dicots_CU1000_onlyEudicots = all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
  filter(max(decimalLatitude) < -23.437)   %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  dplyr::select(genus) %>% 
  distinct(genus) %>% 
  arrange(genus)


# Join all genera
all_genera_CU500 = bind_rows(genera_ntemp_dicots_CU500, 
                             genera_ntemp_monocots_CU500, 
                             genera_stemp_dicots_CU500, 
                             genera_stemp_monocots_CU500, 
                             genera_trop_dicots_CU500, 
                             genera_trop_monocots_CU500)


all_genera_CU500_onlyEudicots = bind_rows(genera_ntemp_dicots_CU500_onlyEudicots, 
                             genera_ntemp_monocots_CU500, 
                             genera_stemp_dicots_CU500_onlyEudicots, 
                             genera_stemp_monocots_CU500, 
                             genera_trop_dicots_CU500_onlyEudicots, 
                             genera_trop_monocots_CU500)

all_genera_CU1000 = bind_rows(genera_ntemp_dicots_CU1000,
                              genera_ntemp_monocots_CU1000,
                              genera_stemp_dicots_CU1000,
                              genera_stemp_monocots_CU1000, 
                              genera_trop_dicots_CU1000, 
                              genera_trop_monocots_CU1000)


all_genera_CU1000_onlyEudicots = bind_rows(genera_ntemp_dicots_CU1000_onlyEudicots,
                              genera_ntemp_monocots_CU1000,
                              genera_stemp_dicots_CU1000_onlyEudicots,
                              genera_stemp_monocots_CU1000, 
                              genera_trop_dicots_CU1000_onlyEudicots, 
                              genera_trop_monocots_CU1000)

# Filter unique genera
all_genera_CU500 = all_genera_CU500 %>% 
  ungroup() %>% 
  distinct(genus)

all_genera_CU500_onlyEudicots = all_genera_CU500_onlyEudicots %>% 
  ungroup() %>% 
  distinct(genus)

all_genera_CU1000 = all_genera_CU1000 %>% 
  ungroup() %>% 
  distinct(genus)

all_genera_CU1000_onlyEudicots = all_genera_CU1000_onlyEudicots %>% 
  ungroup() %>% 
  distinct(genus)


# SISTER SPECIES ANALYSIS ####

#------------------------------------------------------------------------------#
#
# Phylogenetic-based analyses of sister species 
#
# Use the phylogeny provided in V.PhyloMaker to filter for all the genera in our final datasets and then extract species pairs. This tree is the extended version which has the nodes annotated with the node names from the nodes.info.2 table
#
#
#------------------------------------------------------------------------------#

# Read tree with node annotations ####
intree = read.tree("GBOTB.extended.with.assigned.node.name.tre")
#intree = multi2di(intree) #only run if needed; we don't need to run it becuase we are excluding non bifurcating trees. 

# make tibble with node info 
nodes.info.2_tbl = as_tibble(nodes.info.2)

# Define function to get all species pairs ----
get_species_pairs = function(x){
  node = x
  subclade = extract.clade(intree, node)
  if (is.binary(subclade)){
    sp_pairs = extract_sisters(subclade, sis_age = TRUE)
    sp_pairs = as_tibble(sp_pairs) %>%
      mutate_if(is.factor, as.character)
  }
}


# CU 500 #####

# filter tibble of all nodes with genera in our data sets
nodes_genera_CU500 = nodes.info.2_tbl %>%
  filter(genus %in% all_genera_CU500$genus) %>%
  filter(sp.n > 1) %>%
  arrange(genus)
# Resulting object has 1517 genera


# The following lines are commented because the loop to clip subtrees and extract sister species takes a while. We have saved the resulting table as a RDS file which can be read below. Uncomment if you want to run the full analysis.

# Loop through all nodes to extract species pairs CU 500
# sp_pairs_list_CU500 = list()
# for (i in nodes_genera_CU500$rn){
#   tibble_sp_pairs = get_species_pairs(i)
#   sp_pairs_list_CU500[[i]] = tibble_sp_pairs
# }
# 
# # # Bind all rows for all pairs for both tables
# all_species_pairs_CU500 = dplyr::bind_rows(sp_pairs_list_CU500)
# 
# # Fix underscores to compare with other tables
# all_species_pairs_CU500 = all_species_pairs_CU500 %>%
#   mutate(across("sp1", str_replace, "_", " ")) %>%
#   mutate(across("sp2", str_replace, "_", " "))
# 
# 
# # Because we filtered by root nodes, nodes can be shared, therefore subtrees could be pruned multiple times resulting in duplicated records for the pairs. Remove duplicates
# all_species_pairs_CU500 = all_species_pairs_CU500 %>%
#   distinct(pair_age, sp1, sp2)
# # Resulting object has 18844 rows (NOTE: ROWS ARE SPECIES PAIRS)
# 
# saveRDS(all_species_pairs_CU500, "all_species_pairs_CU500.RDS")

all_species_pairs_CU500 = readRDS("all_species_pairs_CU500.RDS")


# Filter species pairs to match species in our datasets (i.e., GBIF records from mountains)
all_species_pairs_CU500_overlap = all_species_pairs_CU500 %>% 
  filter(sp1 %in% all_angiosperms_CU500_allData$species | 
           sp2 %in% all_angiosperms_CU500_allData$species) # and/or filter by pairs
# resulting object has 2012 rows (NOTE: ROWS ARE SPECIES PAIRS)

# Filter species in our dataset that have a match in the list of species pairs for either sp1 and/or sp2
all_angiosperms_CU500_allData_overlap = all_angiosperms_CU500_allData %>% 
  filter(species %in% all_species_pairs_CU500$sp1 | 
           species %in% all_species_pairs_CU500$sp2)
# Resulting object has 2303 rows (NOTE: ROWS ARE SPECIES)

# Check what species from the species pairs are not in our data to look for them in GBIF

# First make the tibbles of pairs a "longer" tibble
all_species_pairs_CU500_overlap_long = pivot_longer(all_species_pairs_CU500_overlap,
                                                    cols = starts_with("sp"),
                                                    names_to = "pair",
                                                    values_to = "species")

# check species not present
species_pairs_CU500_not_inData = all_species_pairs_CU500_overlap_long %>% 
  filter(species %notin% all_angiosperms_CU500_allData_overlap$species)


# Filter sister species ####

# Filter sister species pairs restricted to our dataset. This means species in our dataset that have a match in the list of species pairs for both sp1 and sp2. The outcome of this is that our analysis will be narrowly focused on just sister pairs where each species in the pair:
# - has >=10 records in GBIF
# - has a min elevation of 500 m
# - is not widespread

all_species_pairs_match_all_angiosperms_CU500 = all_species_pairs_CU500 %>% 
  filter(sp1 %in% all_angiosperms_CU500_allData$species & 
           sp2 %in% all_angiosperms_CU500_allData$species)

# Make the list of pairs above longer and join all data
all_species_pairs_match_all_angiosperms_CU500_long_allData = pivot_longer(
  all_species_pairs_match_all_angiosperms_CU500,
  cols = starts_with("sp"),
  names_to = "pair",
  values_to = "species") %>% 
  left_join(all_angiosperms_CU500_allData, 
            by = "species")

# Summarise data: global summary
all_species_pairs_match_all_angiosperms_CU500_long_allData_counts_means_global = all_species_pairs_match_all_angiosperms_CU500_long_allData %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Summarise data by region and taxonomic group
all_species_pairs_match_all_angiosperms_CU500_long_allData_perGroup_perRegion = all_species_pairs_match_all_angiosperms_CU500_long_allData %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Make tables with all specimens per species for the sister species pairs restricted to our dataset. First, make tables with geographic data
all_species_pairs_CU500_fullGeoData = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
  filter(species %in% all_species_pairs_match_all_angiosperms_CU500_long_allData$species) %>% 
  bind_rows(all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>% 
              filter(species %in% all_species_pairs_match_all_angiosperms_CU500_long_allData$species)) 

# Crate summaries per species keeping in the table:
# - elev range
# - min elev
# - max elev

temp_all_species_pairs_CU500_fullGeoData_perSp = all_species_pairs_CU500_fullGeoData %>% 
  summarise(min_elevation = min(elevation),
            max_elevation = max(elevation),
            elevation_range = max(elevation) - min(elevation),
            n_collections = n()) %>%
  dplyr::select(species,
                min_elevation,
                max_elevation,
                elevation_range,
                n_collections)

# Join temp table with table including details about region, taxonomic group, pair, and time
species_pairs_elevation_range_CU500 = left_join(all_species_pairs_match_all_angiosperms_CU500_long_allData,
                                            temp_all_species_pairs_CU500_fullGeoData_perSp,
                                            by = "species") %>% 
  rename(n_collections = n_collections.x,
         elevation_range = elevation_range.x) %>%
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                n_collections,
                min_elevation,
                max_elevation,
                elevation_range)

species_pairs_elevation_range_CU500

# Make table wider to see all pairs per row: USE THIS TABLE FOR ANALYSES
species_pairs_elevation_range_CU500_wide = pivot_wider(
  species_pairs_elevation_range_CU500,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  n_collections, 
                  min_elevation, 
                  max_elevation, 
                  elevation_range), 
  names_glue = "{pair}_{.value}")

species_pairs_elevation_range_CU500_wide

# Estimate elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap = species_pairs_elevation_range_CU500_wide %>% 
  rowwise() %>%
  mutate(elevation_overlap = (min(sp1_max_elevation, 
                                  sp2_max_elevation) - 
                                max(sp1_min_elevation, 
                                    sp2_min_elevation)) / 
           min(sp1_elevation_range,
               sp2_elevation_range)) %>%
  mutate(elevation_overlap = replace(elevation_overlap, 
                                     elevation_overlap < 0, 
                                     0))

species_pairs_elevation_range_CU500_wide_ElevOverlap

# Summary tables ----

# Summary tibble global 
species_pairs_elevation_range_CU500_wide_ElevOverlap_counts_means_global = species_pairs_elevation_range_CU500_wide_ElevOverlap %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))

# Summary tibble per group and per region
species_pairs_elevation_range_CU500_wide_ElevOverlap_counts_means_perGroup_perRegion = species_pairs_elevation_range_CU500_wide_ElevOverlap %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))



# Eudicots ----

# filter tibble of all nodes with genera in our data sets - only Eudicots
nodes_genera_CU500_onlyEudicots = nodes.info.2_tbl %>%
  filter(genus %in% all_genera_CU500_onlyEudicots$genus) %>%
  filter(sp.n > 1) %>%
  arrange(genus)
# resulting object has 1481 genera

# The following lines are commented because the loop to clip subtrees and extract sister species takes a while. We have saved the resulting table as a RDS file which can be read below. Uncomment if you want to run the full analysis.

# #Loop through all nodes to extract species pairs CU 500
# sp_pairs_list_CU500_onlyEudicots = list()
# for (i in nodes_genera_CU500_onlyEudicots$rn){
#   tibble_sp_pairs = get_species_pairs(i)
#   sp_pairs_list_CU500_onlyEudicots[[i]] = tibble_sp_pairs
# }
# # 
# # # Bind all rows for all pairs for both tables
# all_species_pairs_CU500_onlyEudicots = dplyr::bind_rows(sp_pairs_list_CU500_onlyEudicots)
# 
# # Fix underscores to compare with other tables
# all_species_pairs_CU500_onlyEudicots = all_species_pairs_CU500_onlyEudicots %>%
#   mutate(across("sp1", str_replace, "_", " ")) %>%
#   mutate(across("sp2", str_replace, "_", " "))
# 
# # 
# # # Because we filtered by root nodes, nodes can be shared, therefore subtrees could be pruned multiple times resulting in duplicated records for the pairs. Remove duplicates
# all_species_pairs_CU500_onlyEudicots = all_species_pairs_CU500_onlyEudicots %>%
#   distinct(pair_age, sp1, sp2)
# # # Resulting object has 18200 rows (NOTE: ROWS ARE SPECIES PAIRS)
# # 
#  saveRDS(all_species_pairs_CU500_onlyEudicots, "all_species_pairs_CU500_onlyEudicots.RDS")


all_species_pairs_CU500_onlyEudicots = readRDS("all_species_pairs_CU500_onlyEudicots.RDS")

# Filter species pairs to match species in our datasets (i.e., GBIF records from mountains)
all_species_pairs_CU500_overlap_onlyEudicots = all_species_pairs_CU500_onlyEudicots %>% 
  filter(sp1 %in% all_angiosperms_CU500_allData_onlyEudicots$species | 
           sp2 %in% all_angiosperms_CU500_allData_onlyEudicots$species) # and/or filter by pairs
# resulting object has 1994 rows (NOTE: ROWS ARE SPECIES PAIRS)


# Filter species in our dataset that have a match in the list of species pairs for either sp1 and/or sp2
all_angiosperms_CU500_allData_overlap_onlyEudicots = all_angiosperms_CU500_allData_onlyEudicots %>% 
  filter(species %in% all_species_pairs_CU500_onlyEudicots$sp1 | 
           species %in% all_species_pairs_CU500_onlyEudicots$sp2)
# Resulting object has 2284 rows (NOTE: ROWS ARE SPECIES)

# Check what species from the species pairs are not in our data to look for them in GBIF

# First make the tibbles of pairs a "longer" tibble
all_species_pairs_CU500_overlap_long_onlyEudicots = pivot_longer(all_species_pairs_CU500_overlap_onlyEudicots,
                                                    cols = starts_with("sp"),
                                                    names_to = "pair",
                                                    values_to = "species")

# check species not present
species_pairs_CU500_not_inData_onlyEudicots = all_species_pairs_CU500_overlap_long_onlyEudicots %>% 
  filter(species %notin% all_angiosperms_CU500_allData_overlap_onlyEudicots$species)


# Filter sister species pairs ----

# Filter sister species pairs restricted to our dataset. This means species in our dataset that have a match in the list of species pairs for both sp1 and sp2. The outcome of this is that our analysis will be narrowly focused on just sister pairs where each species in the pair:
# - has >=10 records in GBIF
# - has a min elevation of 500 m
# - is not widespread

all_species_pairs_match_all_angiosperms_CU500_onlyEudicots = all_species_pairs_CU500_onlyEudicots %>% 
  filter(sp1 %in% all_angiosperms_CU500_allData_onlyEudicots$species & 
           sp2 %in% all_angiosperms_CU500_allData_onlyEudicots$species)

# Make the list of pairs above longer and join all data
all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots = pivot_longer(
  all_species_pairs_match_all_angiosperms_CU500_onlyEudicots,
  cols = starts_with("sp"),
  names_to = "pair",
  values_to = "species") %>% 
  left_join(all_angiosperms_CU500_allData_onlyEudicots, 
            by = "species")

# Summarise data: global summary
all_species_pairs_match_all_angiosperms_CU500_long_allData_counts_means_global_onlyEudicots = all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Summarise data by region and taxonomic group
all_species_pairs_match_all_angiosperms_CU500_long_allData_perGroup_perRegion_onlyEudicots = all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Make tables with all specimens per species for the sister species pairs restricted to our dataset. First, make tables with geographic data
all_species_pairs_CU500_fullGeoData_onlyEudicots = all_monocots_CU500_nomissID_unique_ptlsFiltered_tenPlus %>%
  filter(species %in% all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots$species) %>% 
  bind_rows(all_dicots_CU500_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
              filter(species %in% all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots$species)) 

# Crate summaries per species keeping in the table:
# - elev range
# - min elev
# - max elev

temp_all_species_pairs_CU500_fullGeoData_perSp_onlyEudicots = all_species_pairs_CU500_fullGeoData_onlyEudicots %>% 
  summarise(min_elevation = min(elevation),
            max_elevation = max(elevation),
            elevation_range = max(elevation) - min(elevation),
            n_collections = n()) %>%
  dplyr::select(species,
                min_elevation,
                max_elevation,
                elevation_range,
                n_collections)

# Join temp table with table including details about region, taxonomic group, pair, and time.

species_pairs_elevation_range_CU500_onlyEudicots = left_join(all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots,
                                                temp_all_species_pairs_CU500_fullGeoData_perSp_onlyEudicots,
                                                by = "species") %>% 
  rename(n_collections = n_collections.x,
         elevation_range = elevation_range.x) %>%
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                n_collections,
                min_elevation,
                max_elevation,
                elevation_range)

species_pairs_elevation_range_CU500_onlyEudicots

# Make table wider to see all pairs per row: USE THIS TABLE FOR ANALYSES
species_pairs_elevation_range_CU500_wide_onlyEudicots = pivot_wider(
  species_pairs_elevation_range_CU500_onlyEudicots,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  n_collections, 
                  min_elevation, 
                  max_elevation, 
                  elevation_range), 
  names_glue = "{pair}_{.value}")

species_pairs_elevation_range_CU500_wide_onlyEudicots

# Estimate elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots = species_pairs_elevation_range_CU500_wide_onlyEudicots %>% 
  rowwise() %>%
  mutate(elevation_overlap = (min(sp1_max_elevation, 
                                  sp2_max_elevation) - 
                                max(sp1_min_elevation, 
                                    sp2_min_elevation)) / 
           min(sp1_elevation_range,
               sp2_elevation_range)) %>%
  mutate(elevation_overlap = replace(elevation_overlap, 
                                     elevation_overlap < 0, 
                                     0))

species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots


 # Summary tables ----
# Summary tibble global 
species_pairs_elevation_range_CU500_wide_ElevOverlap_counts_means_global_onlyEudicots = species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))

# Summary tibble per group and per region
species_pairs_elevation_range_CU500_wide_ElevOverlap_counts_means_perGroup_perRegion_onlyEudicots = species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))


#-------------#
# CU 1000 #####
#-------------#

# filter tibble of all nodes with genera in our data sets
nodes_genera_CU1000 = nodes.info.2_tbl %>%
  filter(genus %in% all_genera_CU1000$genus) %>%
  filter(sp.n > 1) %>%
  arrange(genus)
# Resulting object has 1748 genera


# The following lines are commented because the loop to clip subtrees and extract sister species takes a while. We have saved the resulting table as a RDS file which can be read below. Uncomment if you want to run the full analysis.

# # Loop through all nodes to extract species pairs CU 1000
# sp_pairs_list_CU1000 = list()
# for (i in nodes_genera_CU1000$rn){
#   tibble_sp_pairs = get_species_pairs(i)
#   sp_pairs_list_CU1000[[i]] = tibble_sp_pairs 
# }
# 
# # Bind all rows for all pairs for both tables
# all_species_pairs_CU1000 = dplyr::bind_rows(sp_pairs_list_CU1000)
# 
# # Fix underscores to compare with other tables 
# all_species_pairs_CU1000 = all_species_pairs_CU1000 %>% 
#   mutate(across("sp1", str_replace, "_", " ")) %>% 
#   mutate(across("sp2", str_replace, "_", " "))
# 
# # Because we filtered by root nodes, nodes can be shared, therefore subtrees could be pruned multiple times resulting in duplicated records for the pairs. Remove duplicates
# all_species_pairs_CU1000 = all_species_pairs_CU1000 %>% 
#   distinct(pair_age, sp1, sp2)
# # Resulting object has 19452 rows (NOTE: ROWS ARE SPECIES PAIRS)
# 
# saveRDS(all_species_pairs_CU1000, "all_species_pairs_CU1000.RDS")

all_species_pairs_CU1000 = readRDS("all_species_pairs_CU1000.RDS")


# Filter species pairs to match species in our datasets (i.e., GBIF records from mountains)
all_species_pairs_CU1000_overlap = all_species_pairs_CU1000 %>% 
  filter(sp1 %in% all_angiosperms_CU1000_allData$species | 
           sp2 %in% all_angiosperms_CU1000_allData$species) 
# Resulting has 2774 rows (NOTE: ROWS ARE SPECIES PAIRS)

# Filter species in our dataset that have a match in the list of species pairs for either sp1 and/or sp2
all_angiosperms_CU1000_allData_overlap = all_angiosperms_CU1000_allData %>% 
  filter(species %in% all_species_pairs_CU1000$sp1 | 
           species %in% all_species_pairs_CU1000$sp2)
# Resulting object has 3295 rows (NOTE: ROWS ARE SPECIES))


# Check what species from the species pairs are not in our data to look for them in GBIF

# First make the tibbles of pairs a "longer" tibble
all_species_pairs_CU1000_overlap_long = pivot_longer(all_species_pairs_CU1000_overlap,
                                                    cols=starts_with("sp"),
                                                    names_to = "pair",
                                                    values_to = "species")

# check species not present
species_pairs_CU1000_not_inData = all_species_pairs_CU1000_overlap_long %>% 
  filter(species %notin% all_angiosperms_CU1000_allData_overlap$species)


# Filter sister species pairs ----

# Filter sister species pairs restricted to our dataset. This means species in our dataset that have a match in the list of species pairs for both sp1 and sp2. The outcome of this is that our analysis will be narrowly focused on just sister pairs where each species in the pair:
# - has >=10 records in GBIF
# - has a min elevation of 500 m
# - is not widespread

all_species_pairs_match_all_angiosperms_CU1000 = all_species_pairs_CU1000 %>% 
  filter(sp1 %in% all_angiosperms_CU1000_allData$species & 
           sp2 %in% all_angiosperms_CU1000_allData$species)

# Make the list of pairs above longer and join all data
all_species_pairs_match_all_angiosperms_CU1000_long_allData = pivot_longer(
  all_species_pairs_match_all_angiosperms_CU1000,
  cols = starts_with("sp"),
  names_to = "pair",
  values_to = "species") %>% 
  left_join(all_angiosperms_CU1000_allData, 
            by = "species")

# Summarise data: global summary
all_species_pairs_match_all_angiosperms_CU1000_long_allData_counts_means_global = all_species_pairs_match_all_angiosperms_CU1000_long_allData %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Summarise data by region and taxonomic group
all_species_pairs_match_all_angiosperms_CU1000_long_allData_perGroup_perRegion = all_species_pairs_match_all_angiosperms_CU1000_long_allData %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Make tables with all specimens per species for the sister species pairs restricted to our dataset. First, make tables with geographic data
all_species_pairs_CU1000_fullGeoData = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
  filter(species %in% all_species_pairs_match_all_angiosperms_CU1000_long_allData$species) %>% 
  bind_rows(all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>% 
              filter(species %in% all_species_pairs_match_all_angiosperms_CU1000_long_allData$species)) 


# Crate summaries per species keeping in the table:
# - elev range
# - min elev
# - max elev

temp_all_species_pairs_CU1000_fullGeoData_perSp = all_species_pairs_CU1000_fullGeoData %>% 
  summarise(min_elevation = min(elevation),
         max_elevation = max(elevation),
         elevation_range = max(elevation) - min(elevation),
         n_collections = n()) %>%
  dplyr::select(species,
         min_elevation,
         max_elevation,
         elevation_range,
         n_collections)

# Join temp table with table including details about region, taxonomic group, pair, and time.

species_pairs_elevation_range_CU1000 = left_join(all_species_pairs_match_all_angiosperms_CU1000_long_allData,
                                      temp_all_species_pairs_CU1000_fullGeoData_perSp,
                                      by = "species") %>% 
  rename(n_collections = n_collections.x,
         elevation_range = elevation_range.x) %>%
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                n_collections,
                min_elevation,
                max_elevation,
                elevation_range)

species_pairs_elevation_range_CU1000

# Make table wider to see all pairs per row: USE THIS TABLE FOR ANALYSES
species_pairs_elevation_range_CU1000_wide = pivot_wider(
  species_pairs_elevation_range_CU1000,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  n_collections, 
                  min_elevation, 
                  max_elevation, 
                  elevation_range), 
  names_glue = "{pair}_{.value}")

species_pairs_elevation_range_CU1000_wide

# Estimate elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap = species_pairs_elevation_range_CU1000_wide %>% 
  rowwise() %>%
  mutate(elevation_overlap = (min(sp1_max_elevation, 
                                  sp2_max_elevation) - 
                                max(sp1_min_elevation, 
                                    sp2_min_elevation)) / 
           min(sp1_elevation_range, 
               sp2_elevation_range)) %>%
  mutate(elevation_overlap = replace(elevation_overlap, 
                                     elevation_overlap < 0, 
                                     0))

species_pairs_elevation_range_CU1000_wide_ElevOverlap


# Summary tables ----

# Summary tibble global 
species_pairs_elevation_range_CU1000_wide_ElevOverlap_counts_means_global = species_pairs_elevation_range_CU1000_wide_ElevOverlap %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))

# Summary tibble per group and per region
species_pairs_elevation_range_CU1000_wide_ElevOverlap_counts_means_perGroup_perRegion = species_pairs_elevation_range_CU1000_wide_ElevOverlap %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))


#
# Eudicots ----
#

# Eudicots
nodes_genera_CU1000_onlyEudicots = nodes.info.2_tbl %>%
  filter(genus %in% all_genera_CU1000_onlyEudicots$genus) %>%
  filter(sp.n > 1) %>%
  arrange(genus)
# Resulting object has 1717 genera

# The following lines are commented because the loop to clip subtrees and extract sister species takes a while. We have saved the resulting table as a RDS file which can be read below. Uncomment if you want to run the full analysis.

# # Loop through all nodes to extract species pairs CU 1000
# sp_pairs_list_CU1000_onlyEudicots = list()
# for (i in nodes_genera_CU1000_onlyEudicots$rn){
#   tibble_sp_pairs = get_species_pairs(i)
#   sp_pairs_list_CU1000_onlyEudicots[[i]] = tibble_sp_pairs
# }
# 
# # Bind all rows for all pairs for both tables
# all_species_pairs_CU1000_onlyEudicots = dplyr::bind_rows(sp_pairs_list_CU1000_onlyEudicots)
# 
# # Fix underscores to compare with other tables
# all_species_pairs_CU1000_onlyEudicots = all_species_pairs_CU1000_onlyEudicots %>%
#   mutate(across("sp1", str_replace, "_", " ")) %>%
#   mutate(across("sp2", str_replace, "_", " "))
# 
# # Because we filtered by root nodes, nodes can be shared, therefore subtrees could be pruned multiple times resulting in duplicated records for the pairs. Remove duplicates
# all_species_pairs_CU1000_onlyEudicots = all_species_pairs_CU1000_onlyEudicots %>%
#   distinct(pair_age, sp1, sp2)
# # Resulting object has 18796 rows (NOTE: ROWS ARE SPECIES PAIRS)
# 
# saveRDS(all_species_pairs_CU1000_onlyEudicots, "all_species_pairs_CU1000_onlyEudicots.RDS")

all_species_pairs_CU1000_onlyEudicots = readRDS("all_species_pairs_CU1000_onlyEudicots.RDS")



# Filter species pairs to match species in our datasets (i.e., GBIF records from mountains)
all_species_pairs_CU1000_overlap_onlyEudicots = all_species_pairs_CU1000_onlyEudicots %>% 
  filter(sp1 %in% all_angiosperms_CU1000_allData_onlyEudicots$species | 
           sp2 %in% all_angiosperms_CU1000_allData_onlyEudicots$species) 
# Resulting has 2747 rows (NOTE: ROWS ARE SPECIES PAIRS)

# Filter species in our dataset that have a match in the list of species pairs for either sp1 and/or sp2
all_angiosperms_CU1000_allData_overlap_onlyEudicots = all_angiosperms_CU1000_allData_onlyEudicots %>% 
  filter(species %in% all_species_pairs_CU1000_onlyEudicots$sp1 | 
           species %in% all_species_pairs_CU1000_onlyEudicots$sp2)
# Resulting object has 3266 rows (NOTE: ROWS ARE SPECIES))


# Check what species from the species pairs are not in our data to look for them in GBIF

# First make the tibbles of pairs a "longer" tibble
all_species_pairs_CU1000_overlap_long_onlyEudicots = pivot_longer(all_species_pairs_CU1000_overlap_onlyEudicots,
                                                     cols=starts_with("sp"),
                                                     names_to = "pair",
                                                     values_to = "species")


# check species not present
species_pairs_CU1000_not_inData_onlyEudicots = all_species_pairs_CU1000_overlap_long_onlyEudicots %>% 
  filter(species %notin% all_angiosperms_CU1000_allData_overlap_onlyEudicots$species)



# Filter sister species ----

# Filter sister species pairs restricted to our dataset. This means species in our dataset that have a match in the list of species pairs for both sp1 and sp2. The outcome of this is that our analysis will be narrowly focused on just sister pairs where each species in the pair:
# - has >=10 records in GBIF
# - has a min elevation of 500 m
# - is not widespread

all_species_pairs_match_all_angiosperms_CU1000_onlyEudicots = all_species_pairs_CU1000_onlyEudicots %>% 
  filter(sp1 %in% all_angiosperms_CU1000_allData_onlyEudicots$species & 
           sp2 %in% all_angiosperms_CU1000_allData_onlyEudicots$species)

# Make the list of pairs above longer and join all data
all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots = pivot_longer(
  all_species_pairs_match_all_angiosperms_CU1000_onlyEudicots,
  cols = starts_with("sp"),
  names_to = "pair",
  values_to = "species") %>% 
  left_join(all_angiosperms_CU1000_allData_onlyEudicots, 
            by = "species")

# Summarise data: global summary
all_species_pairs_match_all_angiosperms_CU1000_long_allData_counts_means_global_onlyEudicots = all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots %>% 
  group_by(region) %>% 
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Summarise data by region and taxonomic group
all_species_pairs_match_all_angiosperms_CU1000_long_allData_perGroup_perRegion_onlyEudicots = all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots %>%
  group_by(region, group) %>%
  summarise(n = n(),
            mean_elev = mean(elevation_range),
            mean_temp = mean(temp_range),
            mean_precip = mean(precip_range))

# Make tables with all specimens per species for the sister species pairs restricted to our dataset. First, make tables with geographic data
all_species_pairs_CU1000_fullGeoData_onlyEudicots = all_monocots_CU1000_nomissID_unique_ptlsFiltered_tenPlus %>%
  filter(species %in% all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots$species) %>% 
  bind_rows(all_dicots_CU1000_nomissID_unique_ptlsFiltered_tenPlus_onlyEudicots %>% 
              filter(species %in% all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots$species)) 


# Crate summaries per species keeping in the table:
# - elev range
# - min elev
# - max elev

temp_all_species_pairs_CU1000_fullGeoData_perSp_onlyEudicots = all_species_pairs_CU1000_fullGeoData_onlyEudicots %>% 
  summarise(min_elevation = min(elevation),
            max_elevation = max(elevation),
            elevation_range = max(elevation) - min(elevation),
            n_collections = n()) %>%
  dplyr::select(species,
                min_elevation,
                max_elevation,
                elevation_range,
                n_collections)

# Join temp table with table including details about region, taxonomic group, pair, and time.

species_pairs_elevation_range_CU1000_onlyEudicots = left_join(all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots,
                                                 temp_all_species_pairs_CU1000_fullGeoData_perSp_onlyEudicots,
                                                 by = "species") %>% 
  rename(n_collections = n_collections.x,
         elevation_range = elevation_range.x) %>%
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                n_collections,
                min_elevation,
                max_elevation,
                elevation_range)

species_pairs_elevation_range_CU1000_onlyEudicots

# Make table wider to see all pairs per row: USE THIS TABLE FOR ANALYSES
species_pairs_elevation_range_CU1000_wide_onlyEudicots = pivot_wider(
  species_pairs_elevation_range_CU1000_onlyEudicots,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  n_collections, 
                  min_elevation, 
                  max_elevation, 
                  elevation_range), 
  names_glue = "{pair}_{.value}")

species_pairs_elevation_range_CU1000_wide_onlyEudicots

# Estimate elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots = species_pairs_elevation_range_CU1000_wide_onlyEudicots %>% 
  rowwise() %>%
  mutate(elevation_overlap = (min(sp1_max_elevation, 
                                  sp2_max_elevation) - 
                                max(sp1_min_elevation, 
                                    sp2_min_elevation)) / 
           min(sp1_elevation_range, 
               sp2_elevation_range)) %>%
  mutate(elevation_overlap = replace(elevation_overlap, 
                                     elevation_overlap < 0, 
                                     0))

species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots


# Summary tables ----

# Summary tibble global 
species_pairs_elevation_range_CU1000_wide_ElevOverlap_counts_means_global_onlyEudicots = species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))

# Summary tibble per group and per region
species_pairs_elevation_range_CU1000_wide_ElevOverlap_counts_means_perGroup_perRegion_onlyEudicots = species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_elev_overlap = mean(elevation_overlap),
            sd_elev_overlap = sd(elevation_overlap))


# BioClimatic analysis by species pairs (Prediction 2) #######

#-------------------------------------------------------------------------#
#
# BIOCLIM: Bioclimatic analysis by species pairs
#
# Download version 2.1 climate data for 1970-2000. This version was released in January 2020. We downloaded monthly climate data for minimum and maximum temperature at 30 seconds (~1 km2). Each download is a “zip” file containing 12 GeoTiff (.tif) files, one for each month of the year (January is 1; December is 12). 
#
# Unzip the directory and use each raster to extract values
#
#-------------------------------------------------------------------------#

#*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.#
#*.*.*                                *.*.*#
#*.*.*          DO NOT RUN            *.*.*#
#*.*.*       THE CODE BELOW           *.*.*#
#*.*.*                                *.*.*#
#*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.#


#The code below is *only for reference*. DO NOT RUN. To use the data, load the table already with Bioclim variables at the very end; see below.


# # load Min Temp
# tmin_01 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_01.tif") # Min temp Jan
# tmin_02 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_02.tif") # Min temp Feb
# tmin_03 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_03.tif") # Min temp Mar
# tmin_04 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_04.tif") # Min temp Apr
# tmin_05 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_05.tif") # Min temp May
# tmin_06 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_06.tif") # Min temp Jun
# tmin_07 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_07.tif") # Min temp Jul
# tmin_08 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_08.tif") # Min temp Aug
# tmin_09 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_09.tif") # Min temp Sep
# tmin_10 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_10.tif") # Min temp Oct
# tmin_11 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_11.tif") # Min temp Nov
# tmin_12 = raster("wc2.1_30s_tmin/wc2.1_30s_tmin_12.tif") # Min temp Dec
# 
# # load Max Temp
# tmax_01 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_01.tif") # Max temp Jan
# tmax_02 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_02.tif") # Max temp Feb
# tmax_03 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_03.tif") # Max temp Mar
# tmax_04 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_04.tif") # Max temp Apr
# tmax_05 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_05.tif") # Max temp May
# tmax_06 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_06.tif") # Max temp Jun
# tmax_07 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_07.tif") # Max temp Jul
# tmax_08 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_08.tif") # Max temp Aug
# tmax_09 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_09.tif") # Max temp Sep
# tmax_10 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_10.tif") # Max temp Oct
# tmax_11 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_11.tif") # Min temp Nov
# tmax_12 = raster("wc2.1_30s_tmax/wc2.1_30s_tmax_12.tif") # Min temp Dec
# 
# 
# #Check the raster loaded correctly by plotting
# # plot(tmin_07, # Change to other variables
# #     main = "tmin_07", # Change to other variables
# #     xlab = "Longitude",
# #     ylab = "Latitude",
# #     cex.axis = 1.3,
# #     cex.lab = 1.4,
# #     cex.main = 1.5,
# #     col = rev(heat.colors(10))
# #     )
# 
#
#
# CU 500 #
#
#
#
# Make table only with fewer information

# all_species_pairs_CU500_bioclim = all_species_pairs_CU500_fullGeoData %>%
#   dplyr::select(gbifID,
#                 species,
#                 decimalLatitude,
#                 decimalLongitude) %>%
#   arrange(species)
# 
# # Make dataframes from coordinates to be able to extract values from raster layers
# coords_all_species_pairs_CU500 = data.frame(
#   lon = all_species_pairs_CU500_bioclim$decimalLongitude,
#   lat = all_species_pairs_CU500_bioclim$decimalLatitude)
# coordinates(coords_all_species_pairs_CU500) = c("lon",
#                                                 "lat")
# # Plot to make sure it works
# #map()
# #points(coords_all_species_pairs_CU500, pch=16)
# 
# # Extract values
# tmin_01_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_01,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_02_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_02,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_03_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_03,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_04_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_04,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_05_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_05,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_06_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_06,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_07_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_07,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_08_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_08,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_09_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_09,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_10_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_10,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_11_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_11,
#                                                           y = coords_all_species_pairs_CU500)
# tmin_12_val_all_species_pairs_CU500_fullGeoData = extract(x = tmin_12,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_01_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_01,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_02_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_02,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_03_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_03,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_04_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_04,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_05_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_05,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_06_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_06,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_07_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_07,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_08_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_08,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_09_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_09,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_10_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_10,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_11_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_11,
#                                                           y = coords_all_species_pairs_CU500)
# tmax_12_val_all_species_pairs_CU500_fullGeoData = extract(x = tmax_12,
#                                                           y = coords_all_species_pairs_CU500)
# 
# # Join with tibble
# all_species_pairs_CU500_bioclim = cbind(all_species_pairs_CU500_bioclim,
#                                      tmin_01_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_02_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_03_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_04_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_05_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_06_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_07_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_08_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_09_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_10_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_11_val_all_species_pairs_CU500_fullGeoData,
#                                      tmin_12_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_01_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_02_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_03_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_04_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_05_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_06_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_07_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_08_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_09_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_10_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_11_val_all_species_pairs_CU500_fullGeoData,
#                                      tmax_12_val_all_species_pairs_CU500_fullGeoData
#                                      )
# 
# all_species_pairs_CU500_bioclim = all_species_pairs_CU500_bioclim %>%
#   rename(tmin_01 = ...5,
#          tmin_02 = ...6,
#          tmin_03 = ...7,
#          tmin_04 = ...8,
#          tmin_05 = ...9,
#          tmin_06 = ...10,
#          tmin_07 = ...11,
#          tmin_08 = ...12,
#          tmin_09 = ...13,
#          tmin_10 = ...14,
#          tmin_11 = ...15,
#          tmin_12 = ...16,
#          tmax_01 = ...17,
#          tmax_02 = ...18,
#          tmax_03 = ...19,
#          tmax_04 = ...20,
#          tmax_05 = ...21,
#          tmax_06 = ...22,
#          tmax_07 = ...23,
#          tmax_08 = ...24,
#          tmax_09 = ...25,
#          tmax_10 = ...26,
#          tmax_11 = ...27,
#          tmax_12 = ...28
#          )
# 


# Save tibbles as RDS
#saveRDS(all_species_pairs_CU500_bioclim, "all_species_pairs_CU500_bioclim.RDS")


#
# EUDICOTS #
#


# # # Make table only with fewer information
# # 
# all_species_pairs_CU500_bioclim_onlyEudicots = all_species_pairs_CU500_fullGeoData_onlyEudicots %>%
#    dplyr::select(gbifID,
#                  species, 
#                  decimalLatitude,
#                  decimalLongitude) %>% 
#    arrange(species)
# # 
# # # Make dataframes from coordinates to be able to extract values from raster layers
#  coords_all_species_pairs_CU500_onlyEudicots = data.frame(
#    lon = all_species_pairs_CU500_bioclim_onlyEudicots$decimalLongitude,
#    lat = all_species_pairs_CU500_bioclim_onlyEudicots$decimalLatitude)
#  coordinates(coords_all_species_pairs_CU500_onlyEudicots) = c("lon",
#                                                  "lat")
# # # Plot to make sure it works
# # #map()
# # #points(coords_all_species_pairs_CU500, pch=16)
# # 
# # # Extract values
#  tmin_01_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_01,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_02_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_02,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_03_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_03,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_04_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_04,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_05_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_05,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_06_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_06,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_07_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_07,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_08_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_08,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_09_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_09,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_10_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_10,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_11_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_11,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmin_12_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmin_12,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_01_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_01,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_02_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_02,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_03_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_03,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_04_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_04,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_05_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_05,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_06_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_06,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_07_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_07,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_08_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_08,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_09_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_09,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_10_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_10,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_11_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_11,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  tmax_12_val_all_species_pairs_CU500_fullGeoData_onlyEudicots = extract(x = tmax_12,
#                                                            y = coords_all_species_pairs_CU500_onlyEudicots)
#  
# # # Join with tibble
#  all_species_pairs_CU500_bioclim_onlyEudicots = cbind(all_species_pairs_CU500_bioclim_onlyEudicots,
#                                       tmin_01_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_02_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_03_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_04_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_05_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_06_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_07_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_08_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_09_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_10_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_11_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmin_12_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_01_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_02_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_03_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_04_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_05_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_06_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_07_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_08_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_09_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_10_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_11_val_all_species_pairs_CU500_fullGeoData_onlyEudicots,
#                                       tmax_12_val_all_species_pairs_CU500_fullGeoData_onlyEudicots
#                                       )
# # 
#  all_species_pairs_CU500_bioclim_onlyEudicots = all_species_pairs_CU500_bioclim_onlyEudicots %>%
#    rename(tmin_01 = ...5,
#           tmin_02 = ...6,
#           tmin_03 = ...7,
#           tmin_04 = ...8,
#           tmin_05 = ...9,
#           tmin_06 = ...10,
#           tmin_07 = ...11,
#           tmin_08 = ...12,
#           tmin_09 = ...13,
#           tmin_10 = ...14,
#           tmin_11 = ...15,
#           tmin_12 = ...16,
#           tmax_01 = ...17,
#           tmax_02 = ...18,
#           tmax_03 = ...19,
#           tmax_04 = ...20,
#           tmax_05 = ...21,
#           tmax_06 = ...22,
#           tmax_07 = ...23,
#           tmax_08 = ...24,
#           tmax_09 = ...25,
#           tmax_10 = ...26,
#           tmax_11 = ...27,
#           tmax_12 = ...28
#           )
# # 
# 
# 
# # Save tibbles as RDS
# saveRDS(all_species_pairs_CU500_bioclim_onlyEudicots, "all_species_pairs_CU500_bioclim_onlyEudicots.RDS")




#
# CU 1000 #
#

#
# #
# # # Make table only with fewer information
# #
#  all_species_pairs_CU1000_bioclim = all_species_pairs_CU1000_fullGeoData %>%
#    dplyr::select(gbifID,
#                  species,
#                  decimalLatitude,
#                  decimalLongitude) %>%
#    arrange(species)
# 
# # # Make dataframes from coordinates to be able to extract values from raster layers
#  coords_all_species_pairs_CU1000 = data.frame(
#    lon = all_species_pairs_CU1000_bioclim$decimalLongitude,
#    lat = all_species_pairs_CU1000_bioclim$decimalLatitude)
#  coordinates(coords_all_species_pairs_CU1000) = c("lon",
#                                                  "lat")
# # # Plot to make sure it works
# # #map()
# # #points(coords_all_species_pairs_CU500, pch=16)
# #
# # # Extract values
#  tmin_01_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_01,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_02_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_02,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_03_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_03,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_04_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_04,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_05_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_05,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_06_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_06,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_07_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_07,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_08_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_08,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_09_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_09,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_10_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_10,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_11_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_11,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmin_12_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmin_12,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_01_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_01,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_02_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_02,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_03_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_03,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_04_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_04,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_05_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_05,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_06_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_06,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_07_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_07,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_08_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_08,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_09_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_09,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_10_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_10,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_11_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_11,
#                                                            y = coords_all_species_pairs_CU1000)
#  tmax_12_val_all_species_pairs_CU1000_fullGeoData = extract(x = tmax_12,
#                                                            y = coords_all_species_pairs_CU1000)
# #
# # # Join with tibble
#  all_species_pairs_CU1000_bioclim = cbind(all_species_pairs_CU1000_bioclim,
#                                       tmin_01_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_02_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_03_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_04_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_05_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_06_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_07_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_08_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_09_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_10_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_11_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmin_12_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_01_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_02_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_03_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_04_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_05_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_06_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_07_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_08_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_09_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_10_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_11_val_all_species_pairs_CU1000_fullGeoData,
#                                       tmax_12_val_all_species_pairs_CU1000_fullGeoData
#                                       )
# 
#  all_species_pairs_CU1000_bioclim = all_species_pairs_CU1000_bioclim %>%
#    rename(tmin_01 = ...5,
#           tmin_02 = ...6,
#           tmin_03 = ...7,
#           tmin_04 = ...8,
#           tmin_05 = ...9,
#           tmin_06 = ...10,
#           tmin_07 = ...11,
#           tmin_08 = ...12,
#           tmin_09 = ...13,
#           tmin_10 = ...14,
#           tmin_11 = ...15,
#           tmin_12 = ...16,
#           tmax_01 = ...17,
#           tmax_02 = ...18,
#           tmax_03 = ...19,
#           tmax_04 = ...20,
#           tmax_05 = ...21,
#           tmax_06 = ...22,
#           tmax_07 = ...23,
#           tmax_08 = ...24,
#           tmax_09 = ...25,
#           tmax_10 = ...26,
#           tmax_11 = ...27,
#           tmax_12 = ...28
#           )
# #
# 
# # Save tibbles as RDS
# saveRDS(all_species_pairs_CU1000_bioclim, "all_species_pairs_CU1000_bioclim.RDS")



#
# EUDICOTS #
#

# # # Make table only with fewer information
# # 
# all_species_pairs_CU1000_bioclim_onlyEudicots = all_species_pairs_CU1000_fullGeoData_onlyEudicots %>%
#   dplyr::select(gbifID,
#                 species, 
#                 decimalLatitude,
#                 decimalLongitude) %>% 
#   arrange(species)
# # 
# # # Make dataframes from coordinates to be able to extract values from raster layers
# coords_all_species_pairs_CU1000_onlyEudicots = data.frame(
#   lon = all_species_pairs_CU1000_bioclim_onlyEudicots$decimalLongitude,
#   lat = all_species_pairs_CU1000_bioclim_onlyEudicots$decimalLatitude)
# coordinates(coords_all_species_pairs_CU1000_onlyEudicots) = c("lon",
#                                                              "lat")
# # # Plot to make sure it works
# # #map()
# # #points(coords_all_species_pairs_CU500, pch=16)
# # 
# # # Extract values
# tmin_01_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_01,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_02_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_02,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_03_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_03,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_04_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_04,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_05_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_05,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_06_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_06,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_07_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_07,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_08_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_08,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_09_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_09,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_10_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_10,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_11_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_11,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmin_12_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmin_12,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_01_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_01,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_02_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_02,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_03_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_03,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_04_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_04,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_05_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_05,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_06_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_06,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_07_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_07,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_08_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_08,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_09_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_09,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_10_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_10,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_11_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_11,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# tmax_12_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots = extract(x = tmax_12,
#                                                                        y = coords_all_species_pairs_CU1000_onlyEudicots)
# 
# # # Join with tibble
# all_species_pairs_CU1000_bioclim_onlyEudicots = cbind(all_species_pairs_CU1000_bioclim_onlyEudicots,
#                                                      tmin_01_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_02_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_03_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_04_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_05_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_06_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_07_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_08_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_09_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_10_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_11_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmin_12_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_01_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_02_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_03_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_04_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_05_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_06_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_07_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_08_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_09_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_10_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_11_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots,
#                                                      tmax_12_val_all_species_pairs_CU1000_fullGeoData_onlyEudicots
# )
# # 
# all_species_pairs_CU1000_bioclim_onlyEudicots = all_species_pairs_CU1000_bioclim_onlyEudicots %>%
#   rename(tmin_01 = ...5,
#          tmin_02 = ...6,
#          tmin_03 = ...7,
#          tmin_04 = ...8,
#          tmin_05 = ...9,
#          tmin_06 = ...10,
#          tmin_07 = ...11,
#          tmin_08 = ...12,
#          tmin_09 = ...13,
#          tmin_10 = ...14,
#          tmin_11 = ...15,
#          tmin_12 = ...16,
#          tmax_01 = ...17,
#          tmax_02 = ...18,
#          tmax_03 = ...19,
#          tmax_04 = ...20,
#          tmax_05 = ...21,
#          tmax_06 = ...22,
#          tmax_07 = ...23,
#          tmax_08 = ...24,
#          tmax_09 = ...25,
#          tmax_10 = ...26,
#          tmax_11 = ...27,
#          tmax_12 = ...28
#   )
# # 
# 
# 
# # Save tibbles as RDS
# saveRDS(all_species_pairs_CU1000_bioclim_onlyEudicots, "all_species_pairs_CU1000_bioclim_onlyEudicots.RDS")





#************************************************************************#
#************************************************************************#
#************************************************************************#

#         Read in tibble with BioClim  data for all records. 
#        THESE ARE THE DATA THAT SHOULD BE USED TO ESTIMATE RANGES 
#        IN TEMPERATURE AND PRECIPITAION
#                                                                         
#************************************************************************#
#************************************************************************#
#************************************************************************#

#*****************************************************************#
#                                                                 #
#      YOU HAVE TO *RUN* THE FOLLOWING LINES TO LOAD THE DATA     #
#                                                                 #
#*****************************************************************#

all_species_pairs_CU500_bioclim = readRDS("all_species_pairs_CU500_bioclim.RDS")
all_species_pairs_CU1000_bioclim = readRDS("all_species_pairs_CU1000_bioclim.RDS")

all_species_pairs_CU500_bioclim_onlyEudicots = readRDS("all_species_pairs_CU500_bioclim_onlyEudicots.RDS")
all_species_pairs_CU1000_bioclim_onlyEudicots = readRDS("all_species_pairs_CU1000_bioclim_onlyEudicots.RDS")

# Estimate the ambient temperature range each species experiences each month by calculating the mean minimum (mean min.) and maximum (mean max.) across all of its collection sites. Also estimate range using these means. See Kozak and Wiens 2007 PRSb


#
# CU 500 #####
#


all_species_pairs_CU500_bioclim_perSp = all_species_pairs_CU500_bioclim %>% 
  drop_na() %>% 
  summarise(mean_tmin_01 = mean(tmin_01), 
            mean_tmin_02 = mean(tmin_02),
            mean_tmin_03 = mean(tmin_03),
            mean_tmin_04 = mean(tmin_04),
            mean_tmin_05 = mean(tmin_05),
            mean_tmin_06 = mean(tmin_06),
            mean_tmin_07 = mean(tmin_07),
            mean_tmin_08 = mean(tmin_08),
            mean_tmin_09 = mean(tmin_09),
            mean_tmin_10 = mean(tmin_10),
            mean_tmin_11 = mean(tmin_11),
            mean_tmin_12 = mean(tmin_12),
            mean_tmax_01 = mean(tmax_01),
            mean_tmax_02 = mean(tmax_02),
            mean_tmax_03 = mean(tmax_03),
            mean_tmax_04 = mean(tmax_04),
            mean_tmax_05 = mean(tmax_05),
            mean_tmax_06 = mean(tmax_06),
            mean_tmax_07 = mean(tmax_07),
            mean_tmax_08 = mean(tmax_08),
            mean_tmax_09 = mean(tmax_09),
            mean_tmax_10 = mean(tmax_10),
            mean_tmax_11 = mean(tmax_11),
            mean_tmax_12 = mean(tmax_12),
            range_01 = mean_tmax_01 - mean_tmin_01,
            range_02 = mean_tmax_02 - mean_tmin_02,
            range_03 = mean_tmax_03 - mean_tmin_03,
            range_04 = mean_tmax_04 - mean_tmin_04,
            range_05 = mean_tmax_05 - mean_tmin_05,
            range_06 = mean_tmax_06 - mean_tmin_06,
            range_07 = mean_tmax_07 - mean_tmin_07,
            range_08 = mean_tmax_08 - mean_tmin_08,
            range_09 = mean_tmax_09 - mean_tmin_09,
            range_10 = mean_tmax_10 - mean_tmin_10,
            range_11 = mean_tmax_11 - mean_tmin_11,
            range_12 = mean_tmax_12 - mean_tmin_12,
            n_collections = n())  # count number of collections per species

# Make Table with species pairs and bioclim data

species_pairs_bioclim_CU500 = left_join(
  species_pairs_elevation_range_CU500, 
  all_species_pairs_CU500_bioclim_perSp, 
  by = "species") %>% 
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                mean_tmin_01, 
                mean_tmin_02,
                mean_tmin_03,
                mean_tmin_04,
                mean_tmin_05,
                mean_tmin_06,
                mean_tmin_07,
                mean_tmin_08,
                mean_tmin_09,
                mean_tmin_10,
                mean_tmin_11,
                mean_tmin_12,
                mean_tmax_01,
                mean_tmax_02,
                mean_tmax_03,
                mean_tmax_04,
                mean_tmax_05,
                mean_tmax_06,
                mean_tmax_07,
                mean_tmax_08,
                mean_tmax_09,
                mean_tmax_10,
                mean_tmax_11,
                mean_tmax_12,
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
                range_12,
                n_collections.x) %>%
  rename(n_collections = n_collections.x)
  
species_pairs_bioclim_CU500

# Make table wide to have each pair in one row
species_pairs_bioclim_CU500_wide = pivot_wider(
  species_pairs_bioclim_CU500,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  mean_tmin_01, 
                  mean_tmin_02,
                  mean_tmin_03,
                  mean_tmin_04,
                  mean_tmin_05,
                  mean_tmin_06,
                  mean_tmin_07,
                  mean_tmin_08,
                  mean_tmin_09,
                  mean_tmin_10,
                  mean_tmin_11,
                  mean_tmin_12,
                  mean_tmax_01,
                  mean_tmax_02,
                  mean_tmax_03,
                  mean_tmax_04,
                  mean_tmax_05,
                  mean_tmax_06,
                  mean_tmax_07,
                  mean_tmax_08,
                  mean_tmax_09,
                  mean_tmax_10,
                  mean_tmax_11,
                  mean_tmax_12,
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
                  range_12,
                  n_collections), 
  names_glue = "{pair}_{.value}")

species_pairs_bioclim_CU500_wide


# Estimate temperature range overlap ----

# Table with temperature range overlap

species_pairs_bioclim_CU500_wide_TempOverlap = species_pairs_bioclim_CU500_wide %>% 
  rowwise() %>%
  mutate(temp_overlap_jan = 0.5* ((min(sp1_mean_tmax_01, 
                                      sp2_mean_tmax_01) - 
                                    max(sp1_mean_tmin_01, 
                                        sp2_mean_tmin_01)) / 
                                   (sp1_range_01) +
                                   (min(sp1_mean_tmax_01, 
                                        sp2_mean_tmax_01) - 
                                      max(sp1_mean_tmin_01, 
                                          sp2_mean_tmin_01)) / 
                                   (sp2_range_01)),
         temp_overlap_feb = 0.5* ((min(sp1_mean_tmax_02, 
                                       sp2_mean_tmax_02) - 
                                     max(sp1_mean_tmin_02, 
                                         sp2_mean_tmin_02)) / 
                                    (sp1_range_02) +
                                    (min(sp1_mean_tmax_02, 
                                         sp2_mean_tmax_02) - 
                                       max(sp1_mean_tmin_02, 
                                           sp2_mean_tmin_02)) / 
                                    (sp2_range_02)),
         temp_overlap_mar = 0.5 * ((min(sp1_mean_tmax_03, 
                                        sp2_mean_tmax_03) -
                                      max(sp1_mean_tmin_03, 
                                          sp2_mean_tmin_03)) / 
                                     (sp1_range_03) +
                                     (min(sp1_mean_tmax_03, 
                                          sp2_mean_tmax_03) -
                                        max(sp1_mean_tmin_03, 
                                            sp2_mean_tmin_03)) /
                                     (sp2_range_03)),
         temp_overlap_apr = 0.5 * ((min(sp1_mean_tmax_04, 
                                        sp2_mean_tmax_04) - 
                                      max(sp1_mean_tmin_04, 
                                          sp2_mean_tmin_04)) / 
                                     (sp1_range_04) +
                                     (min(sp1_mean_tmax_04, 
                                          sp2_mean_tmax_04) -
                                        max(sp1_mean_tmin_04, 
                                            sp2_mean_tmin_04)) / 
                                     (sp2_range_04)),
         temp_overlap_may = 0.5 * ((min(sp1_mean_tmax_05, 
                                        sp2_mean_tmax_05) -
                                      max(sp1_mean_tmin_05, 
                                          sp2_mean_tmin_05)) / 
                                     (sp1_range_05) +
                                     (min(sp1_mean_tmax_05,
                                          sp2_mean_tmax_05) -
                                        max(sp1_mean_tmin_05, 
                                            sp2_mean_tmin_05)) / 
                                     (sp2_range_05)),
         temp_overlap_jun = 0.5 * ((min(sp1_mean_tmax_06, 
                                        sp2_mean_tmax_06) - 
                                      max(sp1_mean_tmin_06, 
                                          sp2_mean_tmin_06)) / 
                                     (sp1_range_06) +
                                     (min(sp1_mean_tmax_06, 
                                          sp2_mean_tmax_06) - 
                                        max(sp1_mean_tmin_06, 
                                            sp2_mean_tmin_06)) / 
                                     (sp2_range_06)),
         temp_overlap_jul = 0.5 * ((min(sp1_mean_tmax_07, 
                                        sp2_mean_tmax_07) - 
                                      max(sp1_mean_tmin_07, 
                                          sp2_mean_tmin_07)) / 
                                     (sp1_range_07) +
                                    (min(sp1_mean_tmax_07, 
                                         sp2_mean_tmax_07) -
                                       max(sp1_mean_tmin_07, 
                                           sp2_mean_tmin_07)) / 
                                     (sp2_range_07)),
         temp_overlap_aug = 0.5 * ((min(sp1_mean_tmax_08,
                                        sp2_mean_tmax_08) - 
                                      max(sp1_mean_tmin_08, 
                                          sp2_mean_tmin_08)) / 
                                     (sp1_range_08) +
                                    (min(sp1_mean_tmax_08, 
                                         sp2_mean_tmax_08) - 
                                       max(sp1_mean_tmin_08, 
                                           sp2_mean_tmin_08)) / 
                                     (sp2_range_08)),
         temp_overlap_sep = 0.5 * ((min(sp1_mean_tmax_09, 
                                        sp2_mean_tmax_09) -
                                      max(sp1_mean_tmin_09, 
                                          sp2_mean_tmin_09)) / 
                                     (sp1_range_09) +
                                    (min(sp1_mean_tmax_09, 
                                         sp2_mean_tmax_09) - 
                                       max(sp1_mean_tmin_09, 
                                           sp2_mean_tmin_09)) / 
                                     (sp2_range_09)),
         temp_overlap_oct = 0.5 * ((min(sp1_mean_tmax_10, 
                                        sp2_mean_tmax_10) -
                                      max(sp1_mean_tmin_10, 
                                          sp2_mean_tmin_10)) / 
                                     (sp1_range_10) +
                                    (min(sp1_mean_tmax_10, 
                                         sp2_mean_tmax_10) - 
                                       max(sp1_mean_tmin_10, 
                                           sp2_mean_tmin_10)) / 
                                     (sp2_range_10)),
         temp_overlap_nov = 0.5 * ((min(sp1_mean_tmax_11,
                                        sp2_mean_tmax_11) -
                                      max(sp1_mean_tmin_11, 
                                          sp2_mean_tmin_11)) / 
                                     (sp1_range_11) +
                                    (min(sp1_mean_tmax_11, 
                                         sp2_mean_tmax_11) - 
                                       max(sp1_mean_tmin_11, 
                                           sp2_mean_tmin_11)) / 
                                     (sp2_range_11)),
         temp_overlap_dec = 0.5 * ((min(sp1_mean_tmax_12,
                                        sp2_mean_tmax_12) - 
                                      max(sp1_mean_tmin_12, 
                                          sp2_mean_tmin_12)) / 
                                     (sp1_range_12) +
                                    (min(sp1_mean_tmax_12, 
                                         sp2_mean_tmax_12) - 
                                       max(sp1_mean_tmin_12, 
                                           sp2_mean_tmin_12)) / 
                                     (sp2_range_12))) %>%
  mutate(temp_overlap_jan = replace(temp_overlap_jan, 
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

species_pairs_bioclim_CU500_wide_TempOverlap

# Summary tables

# Summary tibble global 
species_pairs_bioclim_CU500_wide_TempOverlap_counts_means_global = species_pairs_bioclim_CU500_wide_TempOverlap %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))

# Summary tibble per group and per region
species_pairs_bioclim_CU500_wide_TempOverlap_counts_means_perGroup_perRegion = species_pairs_bioclim_CU500_wide_TempOverlap %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))


#          
# EUDICOTS ----
#          


all_species_pairs_CU500_bioclim_perSp_onlyEudicots = all_species_pairs_CU500_bioclim_onlyEudicots %>% 
  drop_na() %>% 
  summarise(mean_tmin_01 = mean(tmin_01), 
            mean_tmin_02 = mean(tmin_02),
            mean_tmin_03 = mean(tmin_03),
            mean_tmin_04 = mean(tmin_04),
            mean_tmin_05 = mean(tmin_05),
            mean_tmin_06 = mean(tmin_06),
            mean_tmin_07 = mean(tmin_07),
            mean_tmin_08 = mean(tmin_08),
            mean_tmin_09 = mean(tmin_09),
            mean_tmin_10 = mean(tmin_10),
            mean_tmin_11 = mean(tmin_11),
            mean_tmin_12 = mean(tmin_12),
            mean_tmax_01 = mean(tmax_01),
            mean_tmax_02 = mean(tmax_02),
            mean_tmax_03 = mean(tmax_03),
            mean_tmax_04 = mean(tmax_04),
            mean_tmax_05 = mean(tmax_05),
            mean_tmax_06 = mean(tmax_06),
            mean_tmax_07 = mean(tmax_07),
            mean_tmax_08 = mean(tmax_08),
            mean_tmax_09 = mean(tmax_09),
            mean_tmax_10 = mean(tmax_10),
            mean_tmax_11 = mean(tmax_11),
            mean_tmax_12 = mean(tmax_12),
            range_01 = mean_tmax_01 - mean_tmin_01,
            range_02 = mean_tmax_02 - mean_tmin_02,
            range_03 = mean_tmax_03 - mean_tmin_03,
            range_04 = mean_tmax_04 - mean_tmin_04,
            range_05 = mean_tmax_05 - mean_tmin_05,
            range_06 = mean_tmax_06 - mean_tmin_06,
            range_07 = mean_tmax_07 - mean_tmin_07,
            range_08 = mean_tmax_08 - mean_tmin_08,
            range_09 = mean_tmax_09 - mean_tmin_09,
            range_10 = mean_tmax_10 - mean_tmin_10,
            range_11 = mean_tmax_11 - mean_tmin_11,
            range_12 = mean_tmax_12 - mean_tmin_12,
            n_collections = n())  # count number of collections per species

# Make Table with species pairs and bioclim data

species_pairs_bioclim_CU500_onlyEudicots = left_join(
  species_pairs_elevation_range_CU500_onlyEudicots, 
  all_species_pairs_CU500_bioclim_perSp_onlyEudicots, 
  by = "species") %>% 
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                mean_tmin_01, 
                mean_tmin_02,
                mean_tmin_03,
                mean_tmin_04,
                mean_tmin_05,
                mean_tmin_06,
                mean_tmin_07,
                mean_tmin_08,
                mean_tmin_09,
                mean_tmin_10,
                mean_tmin_11,
                mean_tmin_12,
                mean_tmax_01,
                mean_tmax_02,
                mean_tmax_03,
                mean_tmax_04,
                mean_tmax_05,
                mean_tmax_06,
                mean_tmax_07,
                mean_tmax_08,
                mean_tmax_09,
                mean_tmax_10,
                mean_tmax_11,
                mean_tmax_12,
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
                range_12,
                n_collections.x) %>%
  rename(n_collections = n_collections.x)

species_pairs_bioclim_CU500_onlyEudicots

# Make table wide to have each per in one row
species_pairs_bioclim_CU500_wide_onlyEudicots = pivot_wider(
  species_pairs_bioclim_CU500_onlyEudicots,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  mean_tmin_01, 
                  mean_tmin_02,
                  mean_tmin_03,
                  mean_tmin_04,
                  mean_tmin_05,
                  mean_tmin_06,
                  mean_tmin_07,
                  mean_tmin_08,
                  mean_tmin_09,
                  mean_tmin_10,
                  mean_tmin_11,
                  mean_tmin_12,
                  mean_tmax_01,
                  mean_tmax_02,
                  mean_tmax_03,
                  mean_tmax_04,
                  mean_tmax_05,
                  mean_tmax_06,
                  mean_tmax_07,
                  mean_tmax_08,
                  mean_tmax_09,
                  mean_tmax_10,
                  mean_tmax_11,
                  mean_tmax_12,
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
                  range_12,
                  n_collections), 
  names_glue = "{pair}_{.value}")

species_pairs_bioclim_CU500_wide_onlyEudicots

# Estimate temperature range overlap ====

# Table with temperature range overlap

species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots = species_pairs_bioclim_CU500_wide_onlyEudicots %>% 
  rowwise() %>%
  mutate(temp_overlap_jan = 0.5* ((min(sp1_mean_tmax_01, 
                                       sp2_mean_tmax_01) - 
                                     max(sp1_mean_tmin_01, 
                                         sp2_mean_tmin_01)) / 
                                    (sp1_range_01) +
                                    (min(sp1_mean_tmax_01, 
                                         sp2_mean_tmax_01) - 
                                       max(sp1_mean_tmin_01, 
                                           sp2_mean_tmin_01)) / 
                                    (sp2_range_01)),
         temp_overlap_feb = 0.5* ((min(sp1_mean_tmax_02, 
                                       sp2_mean_tmax_02) - 
                                     max(sp1_mean_tmin_02, 
                                         sp2_mean_tmin_02)) / 
                                    (sp1_range_02) +
                                    (min(sp1_mean_tmax_02, 
                                         sp2_mean_tmax_02) - 
                                       max(sp1_mean_tmin_02, 
                                           sp2_mean_tmin_02)) / 
                                    (sp2_range_02)),
         temp_overlap_mar = 0.5 * ((min(sp1_mean_tmax_03, 
                                        sp2_mean_tmax_03) -
                                      max(sp1_mean_tmin_03, 
                                          sp2_mean_tmin_03)) / 
                                     (sp1_range_03) +
                                     (min(sp1_mean_tmax_03, 
                                          sp2_mean_tmax_03) -
                                        max(sp1_mean_tmin_03, 
                                            sp2_mean_tmin_03)) /
                                     (sp2_range_03)),
         temp_overlap_apr = 0.5 * ((min(sp1_mean_tmax_04, 
                                        sp2_mean_tmax_04) - 
                                      max(sp1_mean_tmin_04, 
                                          sp2_mean_tmin_04)) / 
                                     (sp1_range_04) +
                                     (min(sp1_mean_tmax_04, 
                                          sp2_mean_tmax_04) -
                                        max(sp1_mean_tmin_04, 
                                            sp2_mean_tmin_04)) / 
                                     (sp2_range_04)),
         temp_overlap_may = 0.5 * ((min(sp1_mean_tmax_05, 
                                        sp2_mean_tmax_05) -
                                      max(sp1_mean_tmin_05, 
                                          sp2_mean_tmin_05)) / 
                                     (sp1_range_05) +
                                     (min(sp1_mean_tmax_05,
                                          sp2_mean_tmax_05) -
                                        max(sp1_mean_tmin_05, 
                                            sp2_mean_tmin_05)) / 
                                     (sp2_range_05)),
         temp_overlap_jun = 0.5 * ((min(sp1_mean_tmax_06, 
                                        sp2_mean_tmax_06) - 
                                      max(sp1_mean_tmin_06, 
                                          sp2_mean_tmin_06)) / 
                                     (sp1_range_06) +
                                     (min(sp1_mean_tmax_06, 
                                          sp2_mean_tmax_06) - 
                                        max(sp1_mean_tmin_06, 
                                            sp2_mean_tmin_06)) / 
                                     (sp2_range_06)),
         temp_overlap_jul = 0.5 * ((min(sp1_mean_tmax_07, 
                                        sp2_mean_tmax_07) - 
                                      max(sp1_mean_tmin_07, 
                                          sp2_mean_tmin_07)) / 
                                     (sp1_range_07) +
                                     (min(sp1_mean_tmax_07, 
                                          sp2_mean_tmax_07) -
                                        max(sp1_mean_tmin_07, 
                                            sp2_mean_tmin_07)) / 
                                     (sp2_range_07)),
         temp_overlap_aug = 0.5 * ((min(sp1_mean_tmax_08,
                                        sp2_mean_tmax_08) - 
                                      max(sp1_mean_tmin_08, 
                                          sp2_mean_tmin_08)) / 
                                     (sp1_range_08) +
                                     (min(sp1_mean_tmax_08, 
                                          sp2_mean_tmax_08) - 
                                        max(sp1_mean_tmin_08, 
                                            sp2_mean_tmin_08)) / 
                                     (sp2_range_08)),
         temp_overlap_sep = 0.5 * ((min(sp1_mean_tmax_09, 
                                        sp2_mean_tmax_09) -
                                      max(sp1_mean_tmin_09, 
                                          sp2_mean_tmin_09)) / 
                                     (sp1_range_09) +
                                     (min(sp1_mean_tmax_09, 
                                          sp2_mean_tmax_09) - 
                                        max(sp1_mean_tmin_09, 
                                            sp2_mean_tmin_09)) / 
                                     (sp2_range_09)),
         temp_overlap_oct = 0.5 * ((min(sp1_mean_tmax_10, 
                                        sp2_mean_tmax_10) -
                                      max(sp1_mean_tmin_10, 
                                          sp2_mean_tmin_10)) / 
                                     (sp1_range_10) +
                                     (min(sp1_mean_tmax_10, 
                                          sp2_mean_tmax_10) - 
                                        max(sp1_mean_tmin_10, 
                                            sp2_mean_tmin_10)) / 
                                     (sp2_range_10)),
         temp_overlap_nov = 0.5 * ((min(sp1_mean_tmax_11,
                                        sp2_mean_tmax_11) -
                                      max(sp1_mean_tmin_11, 
                                          sp2_mean_tmin_11)) / 
                                     (sp1_range_11) +
                                     (min(sp1_mean_tmax_11, 
                                          sp2_mean_tmax_11) - 
                                        max(sp1_mean_tmin_11, 
                                            sp2_mean_tmin_11)) / 
                                     (sp2_range_11)),
         temp_overlap_dec = 0.5 * ((min(sp1_mean_tmax_12,
                                        sp2_mean_tmax_12) - 
                                      max(sp1_mean_tmin_12, 
                                          sp2_mean_tmin_12)) / 
                                     (sp1_range_12) +
                                     (min(sp1_mean_tmax_12, 
                                          sp2_mean_tmax_12) - 
                                        max(sp1_mean_tmin_12, 
                                            sp2_mean_tmin_12)) / 
                                     (sp2_range_12))) %>%
  mutate(temp_overlap_jan = replace(temp_overlap_jan, 
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

species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots

# Summary tables ----

# Summary tibble global 
species_pairs_bioclim_CU500_wide_TempOverlap_counts_means_global_onlyEudicots = species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))

# Summary tibble per group and per region
species_pairs_bioclim_CU500_wide_TempOverlap_counts_means_perGroup_perRegion_onlyEudicots = species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overalp = sd(temperature_overlap))


#
# CU 1000 #####
#


all_species_pairs_CU1000_bioclim_perSp = all_species_pairs_CU1000_bioclim %>% 
  drop_na() %>% 
  summarise(mean_tmin_01 = mean(tmin_01), 
            mean_tmin_02 = mean(tmin_02),
            mean_tmin_03 = mean(tmin_03),
            mean_tmin_04 = mean(tmin_04),
            mean_tmin_05 = mean(tmin_05),
            mean_tmin_06 = mean(tmin_06),
            mean_tmin_07 = mean(tmin_07),
            mean_tmin_08 = mean(tmin_08),
            mean_tmin_09 = mean(tmin_09),
            mean_tmin_10 = mean(tmin_10),
            mean_tmin_11 = mean(tmin_11),
            mean_tmin_12 = mean(tmin_12),
            mean_tmax_01 = mean(tmax_01),
            mean_tmax_02 = mean(tmax_02),
            mean_tmax_03 = mean(tmax_03),
            mean_tmax_04 = mean(tmax_04),
            mean_tmax_05 = mean(tmax_05),
            mean_tmax_06 = mean(tmax_06),
            mean_tmax_07 = mean(tmax_07),
            mean_tmax_08 = mean(tmax_08),
            mean_tmax_09 = mean(tmax_09),
            mean_tmax_10 = mean(tmax_10),
            mean_tmax_11 = mean(tmax_11),
            mean_tmax_12 = mean(tmax_12),
            range_01 = mean_tmax_01 - mean_tmin_01,
            range_02 = mean_tmax_02 - mean_tmin_02,
            range_03 = mean_tmax_03 - mean_tmin_03,
            range_04 = mean_tmax_04 - mean_tmin_04,
            range_05 = mean_tmax_05 - mean_tmin_05,
            range_06 = mean_tmax_06 - mean_tmin_06,
            range_07 = mean_tmax_07 - mean_tmin_07,
            range_08 = mean_tmax_08 - mean_tmin_08,
            range_09 = mean_tmax_09 - mean_tmin_09,
            range_10 = mean_tmax_10 - mean_tmin_10,
            range_11 = mean_tmax_11 - mean_tmin_11,
            range_12 = mean_tmax_12 - mean_tmin_12,
            n_collections = n())  # count number of collections per species

# Make Table with species pairs and bioclim data

species_pairs_bioclim_CU1000 = left_join(
  species_pairs_elevation_range_CU1000, 
  all_species_pairs_CU1000_bioclim_perSp, 
  by = "species") %>% 
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                mean_tmin_01, 
                mean_tmin_02,
                mean_tmin_03,
                mean_tmin_04,
                mean_tmin_05,
                mean_tmin_06,
                mean_tmin_07,
                mean_tmin_08,
                mean_tmin_09,
                mean_tmin_10,
                mean_tmin_11,
                mean_tmin_12,
                mean_tmax_01,
                mean_tmax_02,
                mean_tmax_03,
                mean_tmax_04,
                mean_tmax_05,
                mean_tmax_06,
                mean_tmax_07,
                mean_tmax_08,
                mean_tmax_09,
                mean_tmax_10,
                mean_tmax_11,
                mean_tmax_12,
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
                range_12,
                n_collections.x) %>%
  rename(n_collections = n_collections.x)

species_pairs_bioclim_CU1000

# Make table wide to have each per in one row
species_pairs_bioclim_CU1000_wide = pivot_wider(
  species_pairs_bioclim_CU1000,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  mean_tmin_01, 
                  mean_tmin_02,
                  mean_tmin_03,
                  mean_tmin_04,
                  mean_tmin_05,
                  mean_tmin_06,
                  mean_tmin_07,
                  mean_tmin_08,
                  mean_tmin_09,
                  mean_tmin_10,
                  mean_tmin_11,
                  mean_tmin_12,
                  mean_tmax_01,
                  mean_tmax_02,
                  mean_tmax_03,
                  mean_tmax_04,
                  mean_tmax_05,
                  mean_tmax_06,
                  mean_tmax_07,
                  mean_tmax_08,
                  mean_tmax_09,
                  mean_tmax_10,
                  mean_tmax_11,
                  mean_tmax_12,
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
                  range_12,
                  n_collections), 
  names_glue = "{pair}_{.value}")

species_pairs_bioclim_CU1000_wide

# Estimate temperature range overlap ====

# Table with temperature range overlap

species_pairs_bioclim_CU1000_wide_TempOverlap = species_pairs_bioclim_CU1000_wide %>% 
  rowwise() %>%
  mutate(temp_overlap_jan = 0.5* ((min(sp1_mean_tmax_01, 
                                       sp2_mean_tmax_01) - 
                                     max(sp1_mean_tmin_01, 
                                         sp2_mean_tmin_01)) / 
                                    (sp1_range_01) +
                                    (min(sp1_mean_tmax_01, 
                                         sp2_mean_tmax_01) - 
                                       max(sp1_mean_tmin_01, 
                                           sp2_mean_tmin_01)) / 
                                    (sp2_range_01)),
         temp_overlap_feb = 0.5* ((min(sp1_mean_tmax_02, 
                                       sp2_mean_tmax_02) - 
                                     max(sp1_mean_tmin_02, 
                                         sp2_mean_tmin_02)) / 
                                    (sp1_range_02) +
                                    (min(sp1_mean_tmax_02, 
                                         sp2_mean_tmax_02) - 
                                       max(sp1_mean_tmin_02, 
                                           sp2_mean_tmin_02)) / 
                                    (sp2_range_02)),
         temp_overlap_mar = 0.5 * ((min(sp1_mean_tmax_03, 
                                        sp2_mean_tmax_03) -
                                      max(sp1_mean_tmin_03, 
                                          sp2_mean_tmin_03)) / 
                                     (sp1_range_03) +
                                     (min(sp1_mean_tmax_03, 
                                          sp2_mean_tmax_03) -
                                        max(sp1_mean_tmin_03, 
                                            sp2_mean_tmin_03)) /
                                     (sp2_range_03)),
         temp_overlap_apr = 0.5 * ((min(sp1_mean_tmax_04, 
                                        sp2_mean_tmax_04) - 
                                      max(sp1_mean_tmin_04, 
                                          sp2_mean_tmin_04)) / 
                                     (sp1_range_04) +
                                     (min(sp1_mean_tmax_04, 
                                          sp2_mean_tmax_04) -
                                        max(sp1_mean_tmin_04, 
                                            sp2_mean_tmin_04)) / 
                                     (sp2_range_04)),
         temp_overlap_may = 0.5 * ((min(sp1_mean_tmax_05, 
                                        sp2_mean_tmax_05) -
                                      max(sp1_mean_tmin_05, 
                                          sp2_mean_tmin_05)) / 
                                     (sp1_range_05) +
                                     (min(sp1_mean_tmax_05,
                                          sp2_mean_tmax_05) -
                                        max(sp1_mean_tmin_05, 
                                            sp2_mean_tmin_05)) / 
                                     (sp2_range_05)),
         temp_overlap_jun = 0.5 * ((min(sp1_mean_tmax_06, 
                                        sp2_mean_tmax_06) - 
                                      max(sp1_mean_tmin_06, 
                                          sp2_mean_tmin_06)) / 
                                     (sp1_range_06) +
                                     (min(sp1_mean_tmax_06, 
                                          sp2_mean_tmax_06) - 
                                        max(sp1_mean_tmin_06, 
                                            sp2_mean_tmin_06)) / 
                                     (sp2_range_06)),
         temp_overlap_jul = 0.5 * ((min(sp1_mean_tmax_07, 
                                        sp2_mean_tmax_07) - 
                                      max(sp1_mean_tmin_07, 
                                          sp2_mean_tmin_07)) / 
                                     (sp1_range_07) +
                                     (min(sp1_mean_tmax_07, 
                                          sp2_mean_tmax_07) -
                                        max(sp1_mean_tmin_07, 
                                            sp2_mean_tmin_07)) / 
                                     (sp2_range_07)),
         temp_overlap_aug = 0.5 * ((min(sp1_mean_tmax_08,
                                        sp2_mean_tmax_08) - 
                                      max(sp1_mean_tmin_08, 
                                          sp2_mean_tmin_08)) / 
                                     (sp1_range_08) +
                                     (min(sp1_mean_tmax_08, 
                                          sp2_mean_tmax_08) - 
                                        max(sp1_mean_tmin_08, 
                                            sp2_mean_tmin_08)) / 
                                     (sp2_range_08)),
         temp_overlap_sep = 0.5 * ((min(sp1_mean_tmax_09, 
                                        sp2_mean_tmax_09) -
                                      max(sp1_mean_tmin_09, 
                                          sp2_mean_tmin_09)) / 
                                     (sp1_range_09) +
                                     (min(sp1_mean_tmax_09, 
                                          sp2_mean_tmax_09) - 
                                        max(sp1_mean_tmin_09, 
                                            sp2_mean_tmin_09)) / 
                                     (sp2_range_09)),
         temp_overlap_oct = 0.5 * ((min(sp1_mean_tmax_10, 
                                        sp2_mean_tmax_10) -
                                      max(sp1_mean_tmin_10, 
                                          sp2_mean_tmin_10)) / 
                                     (sp1_range_10) +
                                     (min(sp1_mean_tmax_10, 
                                          sp2_mean_tmax_10) - 
                                        max(sp1_mean_tmin_10, 
                                            sp2_mean_tmin_10)) / 
                                     (sp2_range_10)),
         temp_overlap_nov = 0.5 * ((min(sp1_mean_tmax_11,
                                        sp2_mean_tmax_11) -
                                      max(sp1_mean_tmin_11, 
                                          sp2_mean_tmin_11)) / 
                                     (sp1_range_11) +
                                     (min(sp1_mean_tmax_11, 
                                          sp2_mean_tmax_11) - 
                                        max(sp1_mean_tmin_11, 
                                            sp2_mean_tmin_11)) / 
                                     (sp2_range_11)),
         temp_overlap_dec = 0.5 * ((min(sp1_mean_tmax_12,
                                        sp2_mean_tmax_12) - 
                                      max(sp1_mean_tmin_12, 
                                          sp2_mean_tmin_12)) / 
                                     (sp1_range_12) +
                                     (min(sp1_mean_tmax_12, 
                                          sp2_mean_tmax_12) - 
                                        max(sp1_mean_tmin_12, 
                                            sp2_mean_tmin_12)) / 
                                     (sp2_range_12))) %>%
  mutate(temp_overlap_jan = replace(temp_overlap_jan, 
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

species_pairs_bioclim_CU1000_wide_TempOverlap

# Summary tables ----

# Summary tibble global 
species_pairs_bioclim_CU1000_wide_TempOverlap_counts_means_global = species_pairs_bioclim_CU1000_wide_TempOverlap %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))

# Summary tibble per group and per region
species_pairs_bioclim_CU1000_wide_TempOverlap_counts_means_perGroup_perRegion = species_pairs_bioclim_CU1000_wide_TempOverlap %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))



# EUDICOTS ----
         

all_species_pairs_CU1000_bioclim_perSp_onlyEudicots = all_species_pairs_CU1000_bioclim_onlyEudicots %>% 
  drop_na() %>% 
  summarise(mean_tmin_01 = mean(tmin_01), 
            mean_tmin_02 = mean(tmin_02),
            mean_tmin_03 = mean(tmin_03),
            mean_tmin_04 = mean(tmin_04),
            mean_tmin_05 = mean(tmin_05),
            mean_tmin_06 = mean(tmin_06),
            mean_tmin_07 = mean(tmin_07),
            mean_tmin_08 = mean(tmin_08),
            mean_tmin_09 = mean(tmin_09),
            mean_tmin_10 = mean(tmin_10),
            mean_tmin_11 = mean(tmin_11),
            mean_tmin_12 = mean(tmin_12),
            mean_tmax_01 = mean(tmax_01),
            mean_tmax_02 = mean(tmax_02),
            mean_tmax_03 = mean(tmax_03),
            mean_tmax_04 = mean(tmax_04),
            mean_tmax_05 = mean(tmax_05),
            mean_tmax_06 = mean(tmax_06),
            mean_tmax_07 = mean(tmax_07),
            mean_tmax_08 = mean(tmax_08),
            mean_tmax_09 = mean(tmax_09),
            mean_tmax_10 = mean(tmax_10),
            mean_tmax_11 = mean(tmax_11),
            mean_tmax_12 = mean(tmax_12),
            range_01 = mean_tmax_01 - mean_tmin_01,
            range_02 = mean_tmax_02 - mean_tmin_02,
            range_03 = mean_tmax_03 - mean_tmin_03,
            range_04 = mean_tmax_04 - mean_tmin_04,
            range_05 = mean_tmax_05 - mean_tmin_05,
            range_06 = mean_tmax_06 - mean_tmin_06,
            range_07 = mean_tmax_07 - mean_tmin_07,
            range_08 = mean_tmax_08 - mean_tmin_08,
            range_09 = mean_tmax_09 - mean_tmin_09,
            range_10 = mean_tmax_10 - mean_tmin_10,
            range_11 = mean_tmax_11 - mean_tmin_11,
            range_12 = mean_tmax_12 - mean_tmin_12,
            n_collections = n())  # count number of collections per species

# Make Table with species pairs and bioclim data

species_pairs_bioclim_CU1000_onlyEudicots = left_join(
  species_pairs_elevation_range_CU1000_onlyEudicots, 
  all_species_pairs_CU1000_bioclim_perSp_onlyEudicots, 
  by = "species") %>% 
  dplyr::select(pair_age,
                pair,
                species,
                group,
                region,
                mean_tmin_01, 
                mean_tmin_02,
                mean_tmin_03,
                mean_tmin_04,
                mean_tmin_05,
                mean_tmin_06,
                mean_tmin_07,
                mean_tmin_08,
                mean_tmin_09,
                mean_tmin_10,
                mean_tmin_11,
                mean_tmin_12,
                mean_tmax_01,
                mean_tmax_02,
                mean_tmax_03,
                mean_tmax_04,
                mean_tmax_05,
                mean_tmax_06,
                mean_tmax_07,
                mean_tmax_08,
                mean_tmax_09,
                mean_tmax_10,
                mean_tmax_11,
                mean_tmax_12,
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
                range_12,
                n_collections.x) %>%
  rename(n_collections = n_collections.x)

species_pairs_bioclim_CU1000_onlyEudicots

# Make table wide to have each per in one row
species_pairs_bioclim_CU1000_wide_onlyEudicots = pivot_wider(
  species_pairs_bioclim_CU1000_onlyEudicots,
  names_from = pair,
  values_from = c(species, 
                  group, 
                  region, 
                  mean_tmin_01, 
                  mean_tmin_02,
                  mean_tmin_03,
                  mean_tmin_04,
                  mean_tmin_05,
                  mean_tmin_06,
                  mean_tmin_07,
                  mean_tmin_08,
                  mean_tmin_09,
                  mean_tmin_10,
                  mean_tmin_11,
                  mean_tmin_12,
                  mean_tmax_01,
                  mean_tmax_02,
                  mean_tmax_03,
                  mean_tmax_04,
                  mean_tmax_05,
                  mean_tmax_06,
                  mean_tmax_07,
                  mean_tmax_08,
                  mean_tmax_09,
                  mean_tmax_10,
                  mean_tmax_11,
                  mean_tmax_12,
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
                  range_12,
                  n_collections), 
  names_glue = "{pair}_{.value}")

species_pairs_bioclim_CU1000_wide_onlyEudicots

# Estimate temperature range overlap ====

# Table with temperature range overlap

species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots = species_pairs_bioclim_CU1000_wide_onlyEudicots %>% 
  rowwise() %>%
  mutate(temp_overlap_jan = 0.5* ((min(sp1_mean_tmax_01, 
                                       sp2_mean_tmax_01) - 
                                     max(sp1_mean_tmin_01, 
                                         sp2_mean_tmin_01)) / 
                                    (sp1_range_01) +
                                    (min(sp1_mean_tmax_01, 
                                         sp2_mean_tmax_01) - 
                                       max(sp1_mean_tmin_01, 
                                           sp2_mean_tmin_01)) / 
                                    (sp2_range_01)),
         temp_overlap_feb = 0.5* ((min(sp1_mean_tmax_02, 
                                       sp2_mean_tmax_02) - 
                                     max(sp1_mean_tmin_02, 
                                         sp2_mean_tmin_02)) / 
                                    (sp1_range_02) +
                                    (min(sp1_mean_tmax_02, 
                                         sp2_mean_tmax_02) - 
                                       max(sp1_mean_tmin_02, 
                                           sp2_mean_tmin_02)) / 
                                    (sp2_range_02)),
         temp_overlap_mar = 0.5 * ((min(sp1_mean_tmax_03, 
                                        sp2_mean_tmax_03) -
                                      max(sp1_mean_tmin_03, 
                                          sp2_mean_tmin_03)) / 
                                     (sp1_range_03) +
                                     (min(sp1_mean_tmax_03, 
                                          sp2_mean_tmax_03) -
                                        max(sp1_mean_tmin_03, 
                                            sp2_mean_tmin_03)) /
                                     (sp2_range_03)),
         temp_overlap_apr = 0.5 * ((min(sp1_mean_tmax_04, 
                                        sp2_mean_tmax_04) - 
                                      max(sp1_mean_tmin_04, 
                                          sp2_mean_tmin_04)) / 
                                     (sp1_range_04) +
                                     (min(sp1_mean_tmax_04, 
                                          sp2_mean_tmax_04) -
                                        max(sp1_mean_tmin_04, 
                                            sp2_mean_tmin_04)) / 
                                     (sp2_range_04)),
         temp_overlap_may = 0.5 * ((min(sp1_mean_tmax_05, 
                                        sp2_mean_tmax_05) -
                                      max(sp1_mean_tmin_05, 
                                          sp2_mean_tmin_05)) / 
                                     (sp1_range_05) +
                                     (min(sp1_mean_tmax_05,
                                          sp2_mean_tmax_05) -
                                        max(sp1_mean_tmin_05, 
                                            sp2_mean_tmin_05)) / 
                                     (sp2_range_05)),
         temp_overlap_jun = 0.5 * ((min(sp1_mean_tmax_06, 
                                        sp2_mean_tmax_06) - 
                                      max(sp1_mean_tmin_06, 
                                          sp2_mean_tmin_06)) / 
                                     (sp1_range_06) +
                                     (min(sp1_mean_tmax_06, 
                                          sp2_mean_tmax_06) - 
                                        max(sp1_mean_tmin_06, 
                                            sp2_mean_tmin_06)) / 
                                     (sp2_range_06)),
         temp_overlap_jul = 0.5 * ((min(sp1_mean_tmax_07, 
                                        sp2_mean_tmax_07) - 
                                      max(sp1_mean_tmin_07, 
                                          sp2_mean_tmin_07)) / 
                                     (sp1_range_07) +
                                     (min(sp1_mean_tmax_07, 
                                          sp2_mean_tmax_07) -
                                        max(sp1_mean_tmin_07, 
                                            sp2_mean_tmin_07)) / 
                                     (sp2_range_07)),
         temp_overlap_aug = 0.5 * ((min(sp1_mean_tmax_08,
                                        sp2_mean_tmax_08) - 
                                      max(sp1_mean_tmin_08, 
                                          sp2_mean_tmin_08)) / 
                                     (sp1_range_08) +
                                     (min(sp1_mean_tmax_08, 
                                          sp2_mean_tmax_08) - 
                                        max(sp1_mean_tmin_08, 
                                            sp2_mean_tmin_08)) / 
                                     (sp2_range_08)),
         temp_overlap_sep = 0.5 * ((min(sp1_mean_tmax_09, 
                                        sp2_mean_tmax_09) -
                                      max(sp1_mean_tmin_09, 
                                          sp2_mean_tmin_09)) / 
                                     (sp1_range_09) +
                                     (min(sp1_mean_tmax_09, 
                                          sp2_mean_tmax_09) - 
                                        max(sp1_mean_tmin_09, 
                                            sp2_mean_tmin_09)) / 
                                     (sp2_range_09)),
         temp_overlap_oct = 0.5 * ((min(sp1_mean_tmax_10, 
                                        sp2_mean_tmax_10) -
                                      max(sp1_mean_tmin_10, 
                                          sp2_mean_tmin_10)) / 
                                     (sp1_range_10) +
                                     (min(sp1_mean_tmax_10, 
                                          sp2_mean_tmax_10) - 
                                        max(sp1_mean_tmin_10, 
                                            sp2_mean_tmin_10)) / 
                                     (sp2_range_10)),
         temp_overlap_nov = 0.5 * ((min(sp1_mean_tmax_11,
                                        sp2_mean_tmax_11) -
                                      max(sp1_mean_tmin_11, 
                                          sp2_mean_tmin_11)) / 
                                     (sp1_range_11) +
                                     (min(sp1_mean_tmax_11, 
                                          sp2_mean_tmax_11) - 
                                        max(sp1_mean_tmin_11, 
                                            sp2_mean_tmin_11)) / 
                                     (sp2_range_11)),
         temp_overlap_dec = 0.5 * ((min(sp1_mean_tmax_12,
                                        sp2_mean_tmax_12) - 
                                      max(sp1_mean_tmin_12, 
                                          sp2_mean_tmin_12)) / 
                                     (sp1_range_12) +
                                     (min(sp1_mean_tmax_12, 
                                          sp2_mean_tmax_12) - 
                                        max(sp1_mean_tmin_12, 
                                            sp2_mean_tmin_12)) / 
                                     (sp2_range_12))) %>%
  mutate(temp_overlap_jan = replace(temp_overlap_jan, 
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

species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots

# Summary tables ----

# Summary tibble global 
species_pairs_bioclim_CU1000_wide_TempOverlap_counts_means_global_onlyEudicots = species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots %>% 
  group_by(sp1_region) %>% 
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))

# Summary tibble per group and per region
species_pairs_bioclim_CU1000_wide_TempOverlap_counts_means_perGroup_perRegion_onlyEudicots = species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots %>%
  group_by(sp1_region, sp1_group) %>%
  summarise(n = n(),
            mean_temp_overlap = mean(temperature_overlap),
            sd_temp_overlap = sd(temperature_overlap))




# STATISTICAL ANALYSES ####

            #*^*^*^**^**^*^*^*^*^**^*^*^*^*^*^*^*^*#
            #*^*^*^**^**^*^*^*^*^**^*^*^*^*^*^*^*^*#
            #                                      #
            #         STATISTICAL ANALYSES        #
            #                                      #
            #^*^*^**^**^*^*^*^*^**^*^*^*^*^*^*^*^*^#
            #*^*^*^**^**^*^*^*^*^**^*^*^*^*^*^*^*^*#

# Assumptions of paramatric tests ====

#<><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><>#
#                                      #
# TEST ASSUMPTIONS OF PARAMETRIC TESTS #
#                                      #
#<><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><>#

# Organize tibbles to make testing easier

# Separate tibbles by taxonomic group for analysis of overlap

# CU 500
# Elevation overlap
species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots = species_pairs_elevation_range_CU500_wide_ElevOverlap %>% 
  filter(sp1_group == "Monocots")

species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots = species_pairs_elevation_range_CU500_wide_ElevOverlap %>% 
  filter(sp1_group == "Dicots")

species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots = species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots %>% 
  filter(sp1_group == "Eudicots")

# Temperature overlap
species_pairs_bioclim_CU500_wide_TempOverlap_monocots = species_pairs_bioclim_CU500_wide_TempOverlap %>% 
  filter(sp1_group == "Monocots")

species_pairs_bioclim_CU500_wide_TempOverlap_dicots = species_pairs_bioclim_CU500_wide_TempOverlap %>% 
  filter(sp1_group == "Dicots")

species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots = species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots %>% 
  filter(sp1_group == "Eudicots")

# CU 1000
# Elevation overlap
species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots = species_pairs_elevation_range_CU1000_wide_ElevOverlap %>% 
  filter(sp1_group == "Monocots")

species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots = species_pairs_elevation_range_CU1000_wide_ElevOverlap %>% 
  filter(sp1_group == "Dicots")

species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots = species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots %>% 
  filter(sp1_group == "Eudicots")

# Temperature overlap
species_pairs_bioclim_CU1000_wide_TempOverlap_monocots = species_pairs_bioclim_CU1000_wide_TempOverlap %>% 
  filter(sp1_group == "Monocots")

species_pairs_bioclim_CU1000_wide_TempOverlap_dicots = species_pairs_bioclim_CU1000_wide_TempOverlap %>% 
  filter(sp1_group == "Dicots")

species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots = species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots %>% 
  filter(sp1_group == "Eudicots")



# CU 500 ####

### All angiosperms ####

# Elevation range ----
all_angiosperms_CU500_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_angiosperms_CU500_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  106.25 < 2.2e-16 ***
#       8198                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_angiosperms_CU500_tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_angiosperms_CU500_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  544.51 < 2.2e-16 ***
#       8198                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  6.9646 0.001112 **
#       288                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU500_wide_TempOverlap  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  6.9574 0.00112 **
#       288                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## EUDICOTS ----
# Elevation range ----
all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  101.24 < 2.2e-16 ***
#       8085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  543.51 < 2.2e-16 ***
#       8085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  6.9043 0.001179 **
#       287                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Temperature overlap ----
species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   2  7.2344 0.0008606 ***
#       287                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Testing by taxonomic group ####

# Monocots ----
# Elevation range ----
all_monocots_CU500_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_monocots_CU500_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  49.518 < 2.2e-16 ***
#       1604                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_monocots_CU500_tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_monocots_CU500_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  156.15 < 2.2e-16 ***
#       1604                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  0.6043   0.55
#       55                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU500_wide_TempOverlap_monocots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_monocots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  0.2277 0.7971
#       55                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Dicots ----
# Elevation range ----
all_dicots_CU500_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_dicots_CU500_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  71.421 < 2.2e-16 ***
#       6591                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_dicots_CU500_tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_dicots_CU500_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  397.97 < 2.2e-16 ***
#       6591                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  6.9022 0.001227 **
#       230                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU500_wide_TempOverlap_dicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_dicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  0.2277 0.7971
#       230                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Eudicots ----
# Elevation range ----
all_dicots_CU500_elevRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_dicots_CU500_elevRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  66.849 < 2.2e-16 ***
#       6475                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  397.04 < 2.2e-16 ***
#       6478                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  6.8293 0.001316 **
#       229                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  9.4706 0.0001118 ***
#       232                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# CU 1000 #####


# All angiosperms ----
# Elevation range ----
all_angiosperms_CU1000_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_angiosperms_CU1000_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  127.19 < 2.2e-16 ***
#       12007                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_angiosperms_CU1000_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  416.24 < 2.2e-16 ***
#       12007                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  1.5418  0.215
#       518                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU1000_wide_TempOverlap  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  3.2648 0.03899 *
#       518                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## EUDICOTS ----
# Elevation range ----
all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  120.45 < 2.2e-16 ***
#       11868                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  415.57 < 2.2e-16 ***
#       11868                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  1.5095  0.222
#       516                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Temperature overlap ----
species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   2  3.459 0.03219 *
#       516                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### Testing by taxonomic group ####

# Monocots ----
# Elevation range ----
all_monocots_CU1000_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_monocots_CU1000_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  40.486 < 2.2e-16 ***
#       2201                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_monocots_CU1000_tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_monocots_CU1000_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  117.98 < 2.2e-16 ***
#       2201                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap -----
species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  1.1349 0.3256
#       98                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU1000_wide_TempOverlap_monocots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_monocots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  0.1563 0.8556
#       98                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Dicots ----
# Elevation range -----
all_dicots_CU1000_elevRange_coll_group_region %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_dicots_CU1000_elevRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  95.637 < 2.2e-16 ***
#       9803                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_dicots_CU1000tempRange_precipRange_coll_group_region %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_dicots_CU1000_tempRange_precipRange_coll_group_region)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  306.1 < 2.2e-16 ***
#       9803                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  3.182 0.04251 *
#       417                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap ----
species_pairs_bioclim_CU1000_wide_TempOverlap_dicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_dicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  5.1336 0.006273 **
#       417                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Eudicots ----
# Elevation range ----
all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = elevation_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_range ~ region, all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  88.736 < 2.2e-16 ***
#       9664                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature range ----
all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  ggplot(aes(sample = temp_range)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temp_range ~ region, all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  305.19 < 2.2e-16 ***
#       9664                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Elevation overlap ----
species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots  %>% 
  ggplot(aes(sample = elevation_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(elevation_overlap ~ sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value   Pr(>F)   
# group   2  3.1361 0.04448 *
#       415                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Temperature overlap --0--
species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots  %>% 
  ggplot(aes(sample = temperature_overlap)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(temperature_overlap ~ sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots)

# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)   
# group   2  5.3906 0.004884 **
#       415                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# PREDICTION 1: Differences in ranges ####

#<><><><><><><><><><><><#
#<><><><><><><><><><><><#
#                       #
# DIFFERENCES IN RANGES #
#                       #
#<><><><><><><><><><><><#
#<><><><><><><><><><><><>#


# CU 500 #####

# 1.1 Differences in elevation ranges ----

# * . * . * . * . #
# Elevation range #
# * . * . * . * . #


### All angiosperms ####

# Parametric Tests ----
anova_elevRange_CU500_allAngios = aov(elevation_range ~ region,
                                      data = all_angiosperms_CU500_elevRange_coll_group_region)

res_anova_elevRange_CU500_allAngios = Anova(anova_elevRange_CU500_allAngios, type = 2)
res_anova_elevRange_CU500_allAngios
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region    4.704e+08    2  886.52 < 2.2e-16 ***
#   Residuals 2.175e+09 8198                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_allAngios = TukeyHSD(anova_elevRange_CU500_allAngios)
tukey_anova_elevRange_CU500_allAngios
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_angiosperms_CU500_elevRange_coll_group_region)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -470.46331 -500.96211 -439.9645 0e+00
# Tropics-N. Temperate        78.27848   43.66007  112.8969 4e-07
# Tropics-S. Temperate       548.74179  512.91613  584.5674 0e+00


# Non-parametric tests ----
np_anova_elevRange_CU500_allAngios = kruskal.test(elevation_range ~ region,
             data = all_angiosperms_CU500_elevRange_coll_group_region)
np_anova_elevRange_CU500_allAngios
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2411, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU500_allAngios = pairwise.wilcox.test(all_angiosperms_CU500_elevRange_coll_group_region$elevation_range,
                     all_angiosperms_CU500_elevRange_coll_group_region$region, 
                     p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_allAngios
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_elevRange_coll_group_region$elevation_range and all_angiosperms_CU500_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      1.5e-13      < 2e-16     
# 
# P value adjustment method: bonferroni 


# All angiosperms with EUDICOTS ----

# Parametric Tests ----
anova_elevRange_CU500_allAngios_onlyEudicots = aov(elevation_range ~ region,
                                      data = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots)

res_anova_elevRange_CU500_allAngios_onlyEudicots = Anova(anova_elevRange_CU500_allAngios_onlyEudicots, type = 2)
res_anova_elevRange_CU500_allAngios_onlyEudicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region     460688861    2  861.87 < 2.2e-16 ***
#   Residuals 2160800228 8085                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_allAngios_onlyEudicots = TukeyHSD(anova_elevRange_CU500_allAngios_onlyEudicots)
tukey_anova_elevRange_CU500_allAngios_onlyEudicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -467.42925 -498.20575 -436.6527 0e+00
# Tropics-N. Temperate        83.35617   48.26252  118.4498 1e-07
# Tropics-S. Temperate       550.78541  514.37222  587.1986 0e+00


# Non-parametric tests -----
np_anova_elevRange_CU500_allAngios_onlyEudicots = kruskal.test(elevation_range ~ region,
                                                  data = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots)
np_anova_elevRange_CU500_allAngios_onlyEudicots
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2354.6, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU500_allAngios_onlyEudicots = pairwise.wilcox.test(all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots$elevation_range,
                     all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots$region, 
                     p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_allAngios_onlyEudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots$elevation_range and all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      1.7e-14      < 2e-16     
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots ----

# Parametric Tests ----
anova_elevRange_CU500_monocots = aov(elevation_range ~ region,
                                                   data = all_monocots_CU500_elevRange_coll_group_region)

res_anova_elevRange_CU500_monocots = Anova(anova_elevRange_CU500_monocots, type = 2)
res_anova_elevRange_CU500_monocots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region    100706873    2  159.93 < 2.2e-16 ***
#   Residuals 505004411 1604                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_monocots = TukeyHSD(anova_elevRange_CU500_monocots)
tukey_anova_elevRange_CU500_monocots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_monocots_CU500_elevRange_coll_group_region)
# 
# $region
# diff       lwr        upr     p adj
# S. Temperate-N. Temperate -512.3629 -584.8123 -439.91348 0.0000000
# Tropics-N. Temperate       -17.3186 -110.8076   76.17036 0.9011619
# Tropics-S. Temperate       495.0443  402.1737  587.91494 0.0000000


# Non-parametric tests ----
np_anova_elevRange_CU500_monocots = kruskal.test(elevation_range ~ region,
                                                               data = all_monocots_CU500_elevRange_coll_group_region)
np_anova_elevRange_CU500_monocots 
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 447.67, df = 2, p-value < 2.2e-16
wrst_np_anova_elevRange_CU500_monocots  = 
pairwise.wilcox.test(all_monocots_CU500_elevRange_coll_group_region$elevation_range,
                     all_monocots_CU500_elevRange_coll_group_region$region, 
                     p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_monocots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU500_elevRange_coll_group_region$elevation_range and all_monocots_CU500_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      0.84         <2e-16      
# 
# P value adjustment method: bonferroni 


# Dicots ----

# Parametric Tests ----
anova_elevRange_CU500_dicots = aov(elevation_range ~ region,
                                     data = all_dicots_CU500_elevRange_coll_group_region)

res_anova_elevRange_CU500_dicots = Anova(anova_elevRange_CU500_dicots, type = 2)
res_anova_elevRange_CU500_dicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region     376766050    2  746.82 < 2.2e-16 ***
#   Residuals 1662559502 6591                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_dicots = TukeyHSD(anova_elevRange_CU500_dicots)
tukey_anova_elevRange_CU500_dicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_dicots_CU500_elevRange_coll_group_region)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -464.13038 -497.67615 -430.5846     0
# Tropics-N. Temperate        99.76891   62.82909  136.7087     0
# Tropics-S. Temperate       563.89928  525.15478  602.6438     0


# Non-parametric tests ----
np_anova_elevRange_CU500_dicots = kruskal.test(elevation_range ~ region,
                                                 data = all_dicots_CU500_elevRange_coll_group_region)
np_anova_elevRange_CU500_dicots 
# 
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 1967.4, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU500_dicots = 
pairwise.wilcox.test(all_dicots_CU500_elevRange_coll_group_region$elevation_range,
                     all_dicots_CU500_elevRange_coll_group_region$region, 
                     p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_dicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_elevRange_coll_group_region$elevation_range and all_dicots_CU500_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      9.9e-15      < 2e-16     
# 
# P value adjustment method: bonferroni 


# Eudicots ----

# Parametric Tests ----
anova_elevRange_CU500_eudicots = aov(elevation_range ~ region,
                                   data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots)

res_anova_elevRange_CU500_eudicots = Anova(anova_elevRange_CU500_eudicots, type = 2)
res_anova_elevRange_CU500_eudicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region     366923783    2  720.88 < 2.2e-16 ***
#   Residuals 1648637097 6478                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_eudicots = TukeyHSD(anova_elevRange_CU500_eudicots)
tukey_anova_elevRange_CU500_eudicots 
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -460.3723 -494.31151 -426.4330     0
# Tropics-N. Temperate       105.8740   68.33463  143.4133     0
# Tropics-S. Temperate       566.2463  526.73206  605.7605     0


# Non-parametric tests ----
np_anova_elevRange_CU500_eudicots = kruskal.test(elevation_range ~ region,
                                               data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots)
np_anova_elevRange_CU500_eudicots 

# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 1909.7, df = 2, p-value < 2.2e-16

wrst_anova_elevRange_CU500_eudicots = 
pairwise.wilcox.test(all_dicots_CU500_elevRange_coll_group_region_onlyEudicots$elevation_range,
                     all_dicots_CU500_elevRange_coll_group_region_onlyEudicots$region, 
                     p.adjust.method = "bonf")
wrst_anova_elevRange_CU500_eudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_elevRange_coll_group_region_onlyEudicots$elevation_range and all_dicots_CU500_elevRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      8.8e-16      < 2e-16     
# 
# P value adjustment method: bonferroni 


### Analysis with residuals ####

## Perform all the tests again, but this time use residuals of linear regression (elevation_range ~ n_collections) as the response variable.

# First estimate residuals, test for normality./parametric test, then run tests.

# Estimate residuals ----

###  All Angiosperms ####

# Elevation range ----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU500_elevRange_coll_group_region_res = all_angiosperms_CU500_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU500_elevRange_coll_group_region_res = all_angiosperms_CU500_elevRange_coll_group_region %>% 
  ungroup() %>%
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))
#Plot to check
# all_angiosperms_CU500_elevRange_coll_group_region_res %>%
#   ggplot(aes(x = region,
#              y = res_elevation)) +
#   geom_boxplot() +
#   theme_classic()

# Testing for parametric analyses ----

# Elevation range ----
all_angiosperms_CU500_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_angiosperms_CU500_elevRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  92.811 < 2.2e-16 ***
#       8198                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Only with Eudicots ----

# Elevation range
# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()


# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analyses ----

# Elevation range ----
all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  88.633 < 2.2e-16 ***
#       8085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###  By Taxonomic groups ####

# Monocots ----
# Elevation range ----
#UNCOMMENT FR RESIDUASL PER REGION
# all_monocots_CU500_elevRange_coll_group_region_res = all_monocots_CU500_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall" residuals
all_monocots_CU500_elevRange_coll_group_region_res = all_monocots_CU500_elevRange_coll_group_region %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analyses ----

# Elevation range ----
all_monocots_CU500_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_monocots_CU500_elevRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  48.931 < 2.2e-16 ***
#       1604                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Dicots ----
# Elevation range ----
# UNCOMMNET BELOW FOR RESIDUALS PER REGION
# all_dicots_CU500_elevRange_coll_group_region_res = all_dicots_CU500_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()


# These are the overall residuals
all_dicots_CU500_elevRange_coll_group_region_res = all_dicots_CU500_elevRange_coll_group_region %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analyses  ----

# Elevation range ----
all_dicots_CU500_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_dicots_CU500_elevRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  59.206< 2.2e-16 ***
#       6591                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Eudicots ----
# Elevation range ----
# UNCOMMENT BELOW FOR RESIDUALS PER REGION
# all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the overall residuals
all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analyses ----

# Elevation range ----
all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  55.245 < 2.2e-16 ***
#       6478                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ANOVA and KW with residuals ####

# Parametric Tests ----
anova_elevRange_CU500_allAngios_res = aov(res_elevation ~ region,
                                      data = all_angiosperms_CU500_elevRange_coll_group_region_res)

res_anova_elevRange_CU500_allAngios_res = Anova(anova_elevRange_CU500_allAngios_res, type = 2)
res_anova_elevRange_CU500_allAngios_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     446806250    2   834.9 < 2.2e-16 ***
#   Residuals 2193638182 8198                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tukey_anova_elevRange_CU500_allAngios_res = TukeyHSD(anova_elevRange_CU500_allAngios_res)
tukey_anova_elevRange_CU500_allAngios_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_angiosperms_CU500_elevRange_coll_group_region_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -457.82443 -488.45359 -427.1953 0e+00
# Tropics-N. Temperate        77.79563   43.02924  112.5620 5e-07
# Tropics-S. Temperate       535.62006  499.64126  571.5989 0e+00

# Non-parametric tests ----
np_anova_elevRange_CU500_allAngios_res = kruskal.test(res_elevation ~ region,
                                                  data = all_angiosperms_CU500_elevRange_coll_group_region_res)
np_anova_elevRange_CU500_allAngios_res
# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2293.7, df = 2, p-value = 0.0001381

wrst_np_anova_elevRange_CU500_allAngios_res = pairwise.wilcox.test(all_angiosperms_CU500_elevRange_coll_group_region_res$res_elevation,
                                                                 all_angiosperms_CU500_elevRange_coll_group_region_res$region, 
                                                             p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_allAngios_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_elevRange_coll_group_region_res$res_elevation and all_angiosperms_CU500_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      2e-13        <2e-16      
# 
# P value adjustment method: bonferroni 


# All angiosperms with EUDICOTS ----

# Parametric Tests ----
anova_elevRange_CU500_allAngios_onlyEudicots_res = aov(res_elevation ~ region,
                                                       data =                                       all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res)
res_anova_elevRange_CU500_allAngios_onlyEudicots_res = Anova(anova_elevRange_CU500_allAngios_onlyEudicots_res, type = 2)
res_anova_elevRange_CU500_allAngios_onlyEudicots_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     438384922    2  813.42 < 2.2e-16 ***
#   Residuals 2178653341 8085                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_elevRange_CU500_allAngios_onlyEudicots_res = TukeyHSD(anova_elevRange_CU500_allAngios_onlyEudicots_res)
tukey_elevRange_CU500_allAngios_onlyEudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -455.25454 -486.15793 -424.3512 0e+00
# Tropics-N. Temperate        82.88926   47.65094  118.1276 1e-07
# Tropics-S. Temperate       538.14380  501.58049  574.7071 0e+00

# Non-parametric tests ----
np_anova_elevRange_CU500_allAngios_onlyEudicots_res = kruskal.test(res_elevation ~ region,
                                                                   data = all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res)
np_anova_elevRange_CU500_allAngios_onlyEudicots_res

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2244.2, df = 2, p-value = < 2.2e-16

wrst_np_anova_elevRange_CU500_allAngios_onlyEudicots_res = pairwise.wilcox.test(all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res$res_elevation,
                                                                              all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res$region, 
                                                                              p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_allAngios_onlyEudicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res$res_elevation and all_angiosperms_CU500_elevRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      2.3e-14      < 2e-16     
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots ----

# Parametric Tests ----
anova_elevRange_CU500_monocots_res = aov(res_elevation ~ region,
                                         data = all_monocots_CU500_elevRange_coll_group_region_res)

res_anova_elevRange_CU500_monocots_res = Anova(anova_elevRange_CU500_monocots_res, type = 2)
res_anova_elevRange_CU500_monocots_res

# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     99746724    2  158.12 < 2.2e-16 ***
#   Residuals 505929781 1604                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_monocots_res = TukeyHSD(anova_elevRange_CU500_monocots_res)
tukey_anova_elevRange_CU500_monocots_res

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_monocots_CU500_elevRange_coll_group_region_res)
# 
# $region
# diff       lwr        upr     p adj
# S. Temperate-N. Temperate -510.0131 -582.5288 -437.49729 0.0000000
# Tropics-N. Temperate       -17.5761 -111.1507   75.99847 0.8985341
# Tropics-S. Temperate       492.4370  399.4813  585.39265 0.0000000


# Non-parametric tests ----
np_anova_elevRange_CU500_monocots_res = kruskal.test(res_elevation ~ region,
                                                     data = all_monocots_CU500_elevRange_coll_group_region_res)
np_anova_elevRange_CU500_monocots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 442.48, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU500_monocots_res = pairwise.wilcox.test(all_monocots_CU500_elevRange_coll_group_region_res$res_elevation,
                                                                  all_monocots_CU500_elevRange_coll_group_region_res$region, 
                                                                   p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_monocots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU500_elevRange_coll_group_region_res$res_elevation and all_monocots_CU500_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      0.85         <2e-16      
# 
# P value adjustment method: bonferroni


# Dicots ----

# Parametric Tests ----
anova_elevRange_CU500_dicots_res = aov(res_elevation ~ region,
                                       data = all_dicots_CU500_elevRange_coll_group_region_res)

res_anova_elevRange_CU500_dicots_res = Anova(anova_elevRange_CU500_dicots_res, type = 2)
res_anova_elevRange_CU500_dicots_res

# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     354969373    2  696.71 < 2.2e-16 ***
#   Residuals 1679040083 6591                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_dicots_res = TukeyHSD(anova_elevRange_CU500_dicots_res)
tukey_anova_elevRange_CU500_dicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_dicots_CU500_elevRange_coll_group_region_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -449.23211 -482.94374 -415.5205     0
# Tropics-N. Temperate        99.37328   62.25082  136.4957     0
# Tropics-S. Temperate       548.60539  509.66933  587.5415     0

# Non-parametric tests ----
np_anova_elevRange_CU500_dicots_res = kruskal.test(res_elevation ~ region,
                                                   data = all_dicots_CU500_elevRange_coll_group_region_res)
np_anova_elevRange_CU500_dicots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 1856.7, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU500_dicots_res = 
  pairwise.wilcox.test(all_dicots_CU500_elevRange_coll_group_region_res$res_elevation,
                       all_dicots_CU500_elevRange_coll_group_region_res$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU500_dicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_elevRange_coll_group_region_res$res_elevation and all_dicots_CU500_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      1.3e-14      < 2e-16     
# 
# P value adjustment method: bonferroni 


# Eudicots ----

# Parametric Tests ----
anova_elevRange_CU500_eudicots_res = aov(res_elevation ~ region,
                                         data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res)

res_anova_elevRange_CU500_eudicots_res = Anova(anova_elevRange_CU500_eudicots_res, type = 2)
res_anova_elevRange_CU500_eudicots_res

# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     346349829    2     674 < 2.2e-16 ***
#   Residuals 1664426085 6478                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU500_eudicots_res = TukeyHSD(anova_elevRange_CU500_eudicots_res)
tukey_anova_elevRange_CU500_eudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -445.9547 -480.05605 -411.8533     0
# Tropics-N. Temperate       105.4926   67.77388  143.2113     0
# Tropics-S. Temperate       551.4472  511.74428  591.1502     0

# Non-parametric tests ----
np_anova_elevRange_CU500_eudicots_res = kruskal.test(res_elevation ~ region,
                                                     data = all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res)
np_anova_elevRange_CU500_eudicots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 1805.8,, df = 2, p-value = 0.0005555

wrst_anova_elevRange_CU500_eudicots_res = 
  pairwise.wilcox.test(all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res$res_elevation,
                       all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res$region, 
                       p.adjust.method = "bonf")
wrst_anova_elevRange_CU500_eudicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res$res_elevation and all_dicots_CU500_elevRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      1.2e-15      < 2e-16     
# 
# P value adjustment method: bonferroni Pairwise comparisons using Wilcoxon rank sum test 



# 1.2 Differences in thermal ranges ----

# * . * . * . * . * #
# Temperature range #
# * . * . * . * . * #


### All angiosperms ####

# Parametric Tests ----
anova_tempRange_precipRange_CU500_allAngios = aov(temp_range ~ region,
                                                  data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU500_allAngios = Anova(anova_tempRange_precipRange_CU500_allAngios, type = 2)
res_anova_tempRange_precipRange_CU500_allAngios

# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    292481    2    4884 < 2.2e-16 ***
#   Residuals 245471 8198                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_allAngios = TukeyHSD(anova_tempRange_precipRange_CU500_allAngios)
tukey_anova_tempRange_precipRange_CU500_allAngios

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region)
# 
# $region
# diff         lwr         upr     p adj
# S. Temperate-N. Temperate -12.0111928 -12.3351987 -11.6871870 0.0000000
# Tropics-N. Temperate      -12.2256216 -12.5933925 -11.8578508 0.0000000
# Tropics-S. Temperate       -0.2144288  -0.5950248   0.1661672 0.3836375

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_allAngios = kruskal.test(temp_range ~ region,
                                                              data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU500_allAngios

# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 4841.2, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_allAngios = pairwise.wilcox.test(all_angiosperms_CU500_tempRange_precipRange_coll_group_region$temp_range,
                                                                           all_angiosperms_CU500_tempRange_precipRange_coll_group_region$region, 
                                                                           p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_allAngios

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_tempRange_precipRange_coll_group_region$temp_range and all_angiosperms_CU500_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       1           
# 
# P value adjustment method: bonferroni 

# All angiosperms with EUDICOTS ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_allAngios_onlyEudicots = aov(temp_range ~ region,
                                                               data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)

res_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots = Anova(anova_tempRange_precipRange_CU500_allAngios_onlyEudicots, type = 2)
res_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots

# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    286337    2  4766.1 < 2.2e-16 ***
#   Residuals 242866 8085                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots = TukeyHSD(anova_tempRange_precipRange_CU500_allAngios_onlyEudicots)
tukey_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)
# 
# $region
# diff         lwr         upr     p adj
# S. Temperate-N. Temperate -11.9519805 -12.2782643 -11.6256967 0.0000000
# Tropics-N. Temperate      -12.1563457 -12.5283986 -11.7842928 0.0000000
# Tropics-S. Temperate       -0.2043652  -0.5904075   0.1816771 0.4290318


# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots = kruskal.test(temp_range ~ region,
                                                                           data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)
np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots

# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 4772.5, df = 2, p-value < 2.2e-16
wrst_np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots = pairwise.wilcox.test(all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range,
                                                                                        all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$region, 
                                                                                        p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range and all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       1           
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_monocots = aov(temp_range ~ region,
                                                 data = all_monocots_CU500_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU500_monocots = Anova(anova_tempRange_precipRange_CU500_monocots, type = 2)
res_anova_tempRange_precipRange_CU500_monocots
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region     60874    2  973.25 < 2.2e-16 ***
#   Residuals  50163 1604                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_monocots = TukeyHSD(anova_tempRange_precipRange_CU500_monocots)
tukey_anova_tempRange_precipRange_CU500_monocots

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_monocots_CU500_tempRange_precipRange_coll_group_region)
# 
# $region
# diff       lwr           upr     p adj
# S. Temperate-N. Temperate -12.2610198 -12.98309 -11.538947839 0.0000000
# Tropics-N. Temperate      -13.1828887 -14.11465 -12.251124730 0.0000000
# Tropics-S. Temperate       -0.9218689  -1.84747   0.003732565 0.0512057

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_monocots = kruskal.test(temp_range ~ region,
                                                             data = all_monocots_CU500_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU500_monocots 

# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 938.04, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_monocots  = 
  pairwise.wilcox.test(all_monocots_CU500_tempRange_precipRange_coll_group_region$temp_range,
                       all_monocots_CU500_tempRange_precipRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_monocots

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU500_tempRange_precipRange_coll_group_region$temp_range and all_monocots_CU500_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.0056      
# 
# P value adjustment method: bonferroni 

# Dicots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_dicots = aov(temp_range ~ region,
                                               data = all_dicots_CU500_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU500_dicots = Anova(anova_tempRange_precipRange_CU500_dicots, type = 2)
res_anova_tempRange_precipRange_CU500_dicots

# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    232162    2  3931.5 < 2.2e-16 ***
#   Residuals 194606 6591                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_dicots = TukeyHSD(anova_tempRange_precipRange_CU500_dicots)
tukey_anova_tempRange_precipRange_CU500_dicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_dicots_CU500_tempRange_precipRange_coll_group_region)
# 
# $region
# diff        lwr         upr     p adj
# S. Temperate-N. Temperate -11.99282358 -12.355758 -11.6298895 0.0000000
# Tropics-N. Temperate      -12.01512612 -12.414781 -11.6154717 0.0000000
# Tropics-S. Temperate       -0.02230253  -0.441482   0.3968769 0.9914603

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_dicots = kruskal.test(temp_range ~ region,
                                                           data = all_dicots_CU500_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU500_dicots 

# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 3908.9, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_dicots = 
  pairwise.wilcox.test(all_dicots_CU500_tempRange_precipRange_coll_group_region$temp_range,
                       all_dicots_CU500_tempRange_precipRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_dicots

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_tempRange_precipRange_coll_group_region$temp_range and all_dicots_CU500_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.3         
# 
# P value adjustment method: bonferroni 


# Eudicots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_eudicots = aov(temp_range ~ region,
                                                 data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)

res_anova_tempRange_precipRange_CU500_eudicots = Anova(anova_tempRange_precipRange_CU500_eudicots, type = 2)
res_anova_tempRange_precipRange_CU500_eudicots

# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    226062    2  3812.7 < 2.2e-16 ***
#   Residuals 192047 6478                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_eudicots = TukeyHSD(anova_tempRange_precipRange_CU500_eudicots)
tukey_anova_tempRange_precipRange_CU500_eudicots 

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)
# 
# $region
# diff         lwr         upr     p adj
# S. Temperate-N. Temperate -11.91804258 -12.2843482 -11.5517370 0.0000000
# Tropics-N. Temperate      -11.93151616 -12.3366779 -11.5263545 0.0000000
# Tropics-S. Temperate       -0.01347358  -0.4399498   0.4130026 0.9969805

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_eudicots = kruskal.test(temp_range ~ region,
                                                             data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots)
np_anova_tempRange_precipRange_CU500_eudicots 
# 
# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 3839.1, df = 2, p-value < 2.2e-16

wrst_anova_tempRange_precipRange_CU500_eudicots = 
  pairwise.wilcox.test(all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range,
                       all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$region, 
                       p.adjust.method = "bonf")
wrst_anova_tempRange_precipRange_CU500_eudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range and all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.29        
# 
# P value adjustment method: bonferroni 


### Analyses with residuals ####

## Perform all the tests again, but this time use residuals of linear regression (elevation_range ~ n_collections) as the response variable.

# First estimate residuals, test for normality./parametric test, then run tests.

###  All Angiosperms ####

# Temperature range ----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res = all_angiosperms_CU500_tempRange_precipRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res = all_angiosperms_CU500_tempRange_precipRange_coll_group_region %>% 
  ungroup() %>%
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

#Plot to check
# all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res %>%
#   ggplot(aes(x = region,
#              y = res_temp)) +
#   geom_boxplot() +
#   theme_classic()

# Testing for parametric analysis ----

all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  500.72 < 2.2e-16 ***
#       8198                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Only with Eudicots ----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()


# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

# Testing for parametric analysis ----

all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  499.09 < 2.2e-16 ***
#       8085                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###  By Taxonomic groups ####

# Monocots ----
#UNCOMMENT FR RESIDUASL PER REGION
# all_monocots_CU500_tempRange_precipRange_coll_group_region_res = all_monocots_CU500_tempRange_precipRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall" residuals
all_monocots_CU500_tempRange_precipRange_coll_group_region_res = all_monocots_CU500_tempRange_precipRange_coll_group_region %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

# Testing for parametric analysis ----

all_monocots_CU500_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_monocots_CU500_tempRange_precipRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  149.11 < 2.2e-16 ***
#       1604                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dicots ----
# Temperature range ----

# UNCOMMNET BELOW FOR RESIDUALS PER REGION
# all_dicots_CU500_tempRange_precipRange_coll_group_region_res = all_dicots_CU500_tempRange_precipRangecoll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()


# These are the overall residuals
all_dicots_CU500_tempRange_precipRange_coll_group_region_res = all_dicots_CU500_tempRange_precipRange_coll_group_region %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))


# Testing for parametric analysis ----

all_dicots_CU500_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_dicots_CU500_tempRange_precipRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  362.25 < 2.2e-16 ***
#       6591                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Eudicots ----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION
# all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the overall residuals
all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))


# Testing for parametric analysis ----

all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2   360.8 < 2.2e-16 ***
#       6478                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA and KW with residuals ----

# Parametric Tests ----
# All angiosperms ----
anova_tempRange_precipRange_CU500_allAngios_res = aov(res_temp ~ region,
                                                      data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU500_allAngios_res = Anova(anova_tempRange_precipRange_CU500_allAngios_res, type = 2)
res_anova_tempRange_precipRange_CU500_allAngios_res

# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    286718    2  4703.3 < 2.2e-16 ***
#   Residuals 249879 8198                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_allAngios_res = TukeyHSD(anova_tempRange_precipRange_CU500_allAngios_res)
tukey_anova_tempRange_precipRange_CU500_allAngios_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff         lwr          upr     p adj
# S. Temperate-N. Temperate -11.8022865 -12.1291884 -11.47538447 0.0000000
# Tropics-N. Temperate      -12.2336088 -12.6046670 -11.86255067 0.0000000
# Tropics-S. Temperate       -0.4313224  -0.8153204  -0.04732438 0.0230626

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_allAngios_res = kruskal.test(res_temp ~ region,
                                                                  data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU500_allAngios_res
# 
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 4750.6, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_allAngios_res = pairwise.wilcox.test(all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res$res_temp,
                                                                               all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res$region, 
                                                                               p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_allAngios_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res$res_temp and all_angiosperms_CU500_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.23        
# 
# P value adjustment method: bonferroni 

# All angiosperms with EUDICOTS ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res = aov(res_temp ~ region,
                                                                   data =                                       all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)
res_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res = Anova(anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res, type = 2)
res_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res

# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    280576    2  4586.6 < 2.2e-16 ***
#   Residuals 247291 8085                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_tempRange_precipRange_CU500_allAngios_onlyEudicots_res = TukeyHSD(anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res)
tukey_tempRange_precipRange_CU500_allAngios_onlyEudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff         lwr          upr    p adj
# S. Temperate-N. Temperate -11.741034 -12.0702766 -11.41179088 0.000000
# Tropics-N. Temperate      -12.164441 -12.5398678 -11.78901373 0.000000
# Tropics-S. Temperate       -0.423407  -0.8129503  -0.03386368 0.029228

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res = kruskal.test(res_temp ~ region,
                                                                               data = all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)
np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res

# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 4680.4, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res = pairwise.wilcox.test(all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp,
                                                                                            all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$region, 
                                                                                            p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_allAngios_onlyEudicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp and all_angiosperms_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.26        
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_monocots_res = aov(res_temp ~ region,
                                                     data = all_monocots_CU500_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU500_monocots_res = Anova(anova_tempRange_precipRange_CU500_monocots_res, type = 2)
res_anova_tempRange_precipRange_CU500_monocots_res

# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region     60125    2  948.95 < 2.2e-16 ***
#   Residuals  50815 1604                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
tukey_anova_tempRange_precipRange_CU500_monocots_res = TukeyHSD(anova_tempRange_precipRange_CU500_monocots_res)
tukey_anova_tempRange_precipRange_CU500_monocots_res

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_monocots_CU500_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff        lwr         upr     p adj
# S. Temperate-N. Temperate -12.136688 -12.863434 -11.4099422 0.0000000
# Tropics-N. Temperate      -13.196522 -14.134317 -12.2587262 0.0000000
# Tropics-S. Temperate       -1.059833  -1.991426  -0.1282403 0.0209832

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_monocots_res = kruskal.test(res_temp ~ region,
                                                                 data = all_monocots_CU500_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU500_monocots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 926.43, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_monocots_res = pairwise.wilcox.test(all_monocots_CU500_tempRange_precipRange_coll_group_region_res$res_temp,
                                                                              all_monocots_CU500_tempRange_precipRange_coll_group_region_res$region, 
                                                                              p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_monocots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU500_tempRange_precipRange_coll_group_region_res$res_temp and all_monocots_CU500_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.0013      
# 
# P value adjustment method: bonferroni 


# Dicots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_dicots_res = aov(res_temp ~ region,
                                                   data = all_dicots_CU500_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU500_dicots_res = Anova(anova_tempRange_precipRange_CU500_dicots_res, type = 2)
res_anova_tempRange_precipRange_CU500_dicots_res

# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    227191    2  3775.4 < 2.2e-16 ***
#   Residuals 198310 6591                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_dicots_res = TukeyHSD(anova_tempRange_precipRange_CU500_dicots_res)
tukey_anova_tempRange_precipRange_CU500_dicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_dicots_CU500_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff         lwr         upr     p adj
# S. Temperate-N. Temperate -11.7627723 -12.1291434 -11.3964012 0.0000000
# Tropics-N. Temperate      -12.0212407 -12.4246799 -11.6178014 0.0000000
# Tropics-S. Temperate       -0.2584684  -0.6816176   0.1646808 0.3245047

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_dicots_res = kruskal.test(res_temp ~ region,
                                                               data = all_dicots_CU500_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU500_dicots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 3829.7, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU500_dicots_res = 
  pairwise.wilcox.test(all_dicots_CU500_tempRange_precipRange_coll_group_region_res$res_temp,
                       all_dicots_CU500_tempRange_precipRange_coll_group_region_res$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU500_dicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_tempRange_precipRange_coll_group_region_res$res_temp and all_dicots_CU500_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       1           
# 
# P value adjustment method: bonferroni 

# Eudicots ----

# Parametric Tests ----
anova_tempRange_precipRange_CU500_eudicots_res = aov(res_temp ~ region,
                                                     data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)

res_anova_tempRange_precipRange_CU500_eudicots_res = Anova(anova_tempRange_precipRange_CU500_eudicots_res, type = 2)
res_anova_tempRange_precipRange_CU500_eudicots_res

# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    221097    2  3658.2 < 2.2e-16 ***
#   Residuals 195761 6478                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU500_eudicots_res = TukeyHSD(anova_tempRange_precipRange_CU500_eudicots_res)
tukey_anova_tempRange_precipRange_CU500_eudicots_res

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff         lwr         upr     p adj
# S. Temperate-N. Temperate -11.6849416 -12.0547722 -11.3151109 0.0000000
# Tropics-N. Temperate      -11.9376870 -12.3467477 -11.5286263 0.0000000
# Tropics-S. Temperate       -0.2527454  -0.6833257   0.1778349 0.3535642

# Non-parametric tests ----
np_anova_tempRange_precipRange_CU500_eudicots_res = kruskal.test(res_temp ~ region,
                                                                 data = all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res)
np_anova_tempRange_precipRange_CU500_eudicots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 3758.6, df = 2, p-value < 2.2e-16

wrst_anova_tempRange_precipRange_CU500_eudicots_res = 
  pairwise.wilcox.test(all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp,
                       all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$region, 
                       p.adjust.method = "bonf")
wrst_anova_tempRange_precipRange_CU500_eudicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp and all_dicots_CU500_tempRange_precipRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       1           
# 
# P value adjustment method: bonferroni 



#-------------#
# CU 1000 #####
#-------------#

# 1.1 Differences in elevation range ----

# * . * . * . * . #
# Elevation range #
# * . * . * . * . #


### All angiosperms ####

# Parametric Tests ----
anova_elevRange_CU1000_allAngios = aov(elevation_range ~ region,
                                       data = all_angiosperms_CU1000_elevRange_coll_group_region)

res_anova_elevRange_CU1000_allAngios = Anova(anova_elevRange_CU1000_allAngios, type = 2)
res_anova_elevRange_CU1000_allAngios
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq    Df F value    Pr(>F)    
# region     527564898     2  980.35 < 2.2e-16 ***
#   Residuals 3230718599 12007                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
tukey_anova_elevRange_CU1000_allAngios = TukeyHSD(anova_elevRange_CU1000_allAngios)
tukey_anova_elevRange_CU1000_allAngios
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_angiosperms_CU1000_elevRange_coll_group_region)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -457.98264 -484.76792 -431.1974     0
# Tropics-N. Temperate        93.83561   63.98493  123.6863     0
# Tropics-S. Temperate       551.81824  517.48242  586.1541     0

# Non-parametric tests ----
np_anova_elevRange_CU1000_allAngios = kruskal.test(elevation_range ~ region,
                                                   data = all_angiosperms_CU1000_elevRange_coll_group_region)
np_anova_elevRange_CU1000_allAngios
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2739.4, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_allAngios = pairwise.wilcox.test(all_angiosperms_CU1000_elevRange_coll_group_region$elevation_range,
                                                                all_angiosperms_CU1000_elevRange_coll_group_region$region, 
                                                                p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_allAngios
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_elevRange_coll_group_region$elevation_range and all_angiosperms_CU1000_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      1.5e-13      < 2e-16     
# 
# P value adjustment method: bonferroni 

# All angiosperms with EUDICOTS ----

# Parametric Tests ----
anova_elevRange_CU1000_allAngios_onlyEudicots = aov(elevation_range ~ region,
                                                    data = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots)

res_anova_elevRange_CU1000_allAngios_onlyEudicots = Anova(anova_elevRange_CU1000_allAngios_onlyEudicots, type = 2)
res_anova_elevRange_CU1000_allAngios_onlyEudicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq    Df F value    Pr(>F)    
# region     514772241     2  952.41 < 2.2e-16 ***
#   Residuals 3207279835 11868                      
# ---
#   # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_allAngios_onlyEudicots = TukeyHSD(anova_elevRange_CU1000_allAngios_onlyEudicots)
tukey_anova_elevRange_CU1000_allAngios_onlyEudicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -455.02153 -482.04434 -427.9987     0
# Tropics-N. Temperate        97.74961   67.43234  128.0669     0
# Tropics-S. Temperate       552.77114  517.88316  587.6591     0


# Non-parametric tests ----
np_anova_elevRange_CU1000_allAngios_onlyEudicots = kruskal.test(elevation_range ~ region,
                                                                data = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots)
np_anova_elevRange_CU1000_allAngios_onlyEudicots
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2676, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_allAngios_onlyEudicots = pairwise.wilcox.test(all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots$elevation_range,
                                                                             all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots$region, 
                                                                             p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_allAngios_onlyEudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots$elevation_range and all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots----

# Parametric Tests----
anova_elevRange_CU1000_monocots = aov(elevation_range ~ region,
                                      data = all_monocots_CU1000_elevRange_coll_group_region)

res_anova_elevRange_CU1000_monocots = Anova(anova_elevRange_CU1000_monocots, type = 2)
res_anova_elevRange_CU1000_monocots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region    124625119    2  200.47 < 2.2e-16 ***
#   Residuals 684142414 2201                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_monocots = TukeyHSD(anova_elevRange_CU1000_monocots)
tukey_anova_elevRange_CU1000_monocots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_monocots_CU1000_elevRange_coll_group_region)
# 
# $region
# diff       lwr        upr     p adj
# S. Temperate-N. Temperate -516.79675 -579.2682 -454.32528 0.0000000
# Tropics-N. Temperate       -24.07201 -106.2229   58.07886 0.7709831
# Tropics-S. Temperate       492.72474  404.5592  580.89031 0.0000000


# Non-parametric tests----
np_anova_elevRange_CU1000_monocots = kruskal.test(elevation_range ~ region,
                                                  data = all_monocots_CU1000_elevRange_coll_group_region)
np_anova_elevRange_CU1000_monocots 
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 578.28, df = 2, p-value < 2.2e-16
wrst_np_anova_elevRange_CU1000_monocots  = 
  pairwise.wilcox.test(all_monocots_CU1000_elevRange_coll_group_region$elevation_range,
                       all_monocots_CU1000_elevRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_monocots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU1000_elevRange_coll_group_region$elevation_range and all_monocots_CU1000_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      0.93         <2e-16      
# 
# P value adjustment method: bonferroni 


# Dicots----

# Parametric Tests----
anova_elevRange_CU1000_dicots = aov(elevation_range ~ region,
                                    data = all_dicots_CU1000_elevRange_coll_group_region)

res_anova_elevRange_CU1000_dicots = Anova(anova_elevRange_CU1000_dicots, type = 2)
res_anova_elevRange_CU1000_dicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region     418294397    2  811.18 < 2.2e-16 ***
#   Residuals 2527527155 9803                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_dicots = TukeyHSD(anova_elevRange_CU1000_dicots)
tukey_anova_elevRange_CU1000_dicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_dicots_CU1000_elevRange_coll_group_region)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -450.6619 -480.31384 -421.0099     0
# Tropics-N. Temperate       118.0825   86.32112  149.8438     0
# Tropics-S. Temperate       568.7444  531.45492  606.0338     0

# Non-parametric tests----
np_anova_elevRange_CU1000_dicots = kruskal.test(elevation_range ~ region,
                                                data = all_dicots_CU1000_elevRange_coll_group_region)
np_anova_elevRange_CU1000_dicots 
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2197.9, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_dicots = 
  pairwise.wilcox.test(all_dicots_CU1000_elevRange_coll_group_region$elevation_range,
                       all_dicots_CU1000_elevRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_dicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_elevRange_coll_group_region$elevation_range and all_dicots_CU1000_elevRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


# Eudicots----

# Parametric Tests----
anova_elevRange_CU1000_eudicots = aov(elevation_range ~ region,
                                      data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots)

res_anova_elevRange_CU1000_eudicots = Anova(anova_elevRange_CU1000_eudicots, type = 2)
res_anova_elevRange_CU1000_eudicots
# Anova Table (Type II tests)
# 
# Response: elevation_range
# Sum Sq   Df F value    Pr(>F)    
# region     405539768    2  782.43 < 2.2e-16 ***
#   Residuals 2504474601 9664                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_eudicots = TukeyHSD(anova_elevRange_CU1000_eudicots)
tukey_anova_elevRange_CU1000_eudicots 
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = elevation_range ~ region, data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -447.0342 -477.02734 -417.0410     0
# Tropics-N. Temperate       122.8040   90.46094  155.1471     0
# Tropics-S. Temperate       569.8382  531.83115  607.8452     0


# Non-parametric tests----
np_anova_elevRange_CU1000_eudicots = kruskal.test(elevation_range ~ region,
                                                  data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots)
np_anova_elevRange_CU1000_eudicots 
# 
# Kruskal-Wallis rank sum test
# 
# data:  elevation_range by region
# Kruskal-Wallis chi-squared = 2133.3, df = 2, p-value < 2.2e-16

wrst_anova_elevRange_CU1000_eudicots = 
  pairwise.wilcox.test(all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots$elevation_range,
                       all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots$region, 
                       p.adjust.method = "bonf")
wrst_anova_elevRange_CU1000_eudicots

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots$elevation_range and all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


### Analysis with residuals ####

## Perform all the tests again, but this time use residuals of linear regression (elevation_range ~ n_collections) as the response variable.

# First estimate residuals, test for normality./parametric test, then run tests.

# Estimate residuals ----

###  All Angiosperms ####

# Elevation range ----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU1000_elevRange_coll_group_region_res = all_angiosperms_CU1000_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU1000_elevRange_coll_group_region_res = all_angiosperms_CU1000_elevRange_coll_group_region %>% 
  ungroup() %>%
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))
#Plot to check
# all_angiosperms_CU1000_elevRange_coll_group_region_res %>%
#   ggplot(aes(x = region,
#              y = res_elevation)) +
#   geom_boxplot() +
#   theme_classic()

# Testing for parametric analysis ----

all_angiosperms_CU1000_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_angiosperms_CU1000_elevRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
# group     2  126.22 < 2.2e-16 ***
#       12007                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Only with Eudicots----

# Elevation range
# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()


# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))

# Testing for parametric analysis ----

all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
# group     2  120.56 < 2.2e-16 ***
#       11868                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###  By Taxonomic groups ####

# Monocots ----

#UNCOMMENT FOR RESIDUALS PER REGION
# all_monocots_CU1000_elevRange_coll_group_region_res = all_monocots_CU1000_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall" residuals
all_monocots_CU1000_elevRange_coll_group_region_res = all_monocots_CU1000_elevRange_coll_group_region %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analysis ----

all_monocots_CU1000_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_monocots_CU1000_elevRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  42.096 < 2.2e-16 ***
#       2201                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Dicots----

# UNCOMMNET BELOW FOR RESIDUALS PER REGION
# all_dicots_CU1000_elevRange_coll_group_region_res = all_dicots_CU1000_elevRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()


# These are the overall residuals
all_dicots_CU1000_elevRange_coll_group_region_res = all_dicots_CU1000_elevRange_coll_group_region %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# Testing for parametric analysis ----

all_dicots_CU1000_elevRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_dicots_CU1000_elevRange_coll_group_region_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2   91.91 < 2.2e-16 ***
#       9803                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Eudicots----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION
# all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(elevation_range ~ n_collections))) %>% 
#   ungroup()

# These are the overall residuals
all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots %>% 
  mutate(res_elevation = residuals(lm(elevation_range ~ n_collections)))


# ATesting for parametric analysis ----

all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_elevation)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_elevation ~ region, all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res)

# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  86.157 < 2.2e-16 ***
#       9664                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ANOVA and KW with residuals----

# All angiosperms ----

# Parametric Tests----
anova_elevRange_CU1000_allAngios_res = aov(res_elevation ~ region,
                                           data = all_angiosperms_CU1000_elevRange_coll_group_region_res)

res_anova_elevRange_CU1000_allAngios_res = Anova(anova_elevRange_CU1000_allAngios_res, type = 2)
res_anova_elevRange_CU1000_allAngios_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq    Df F value    Pr(>F)    
# region     525690277     2  976.31 < 2.2e-16 ***
#   Residuals 3232560674 12007                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tukey_anova_elevRange_CU1000_allAngios_res = TukeyHSD(anova_elevRange_CU1000_allAngios_res)
tukey_anova_elevRange_CU1000_allAngios_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_angiosperms_CU1000_elevRange_coll_group_region_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -457.18644 -483.97935 -430.3935     0
# Tropics-N. Temperate        93.62167   63.76248  123.4809     0
# Tropics-S. Temperate       550.80811  516.46250  585.1537     0

# Non-parametric tests----
np_anova_elevRange_CU1000_allAngios_res = kruskal.test(res_elevation ~ region,
                                                       data = all_angiosperms_CU1000_elevRange_coll_group_region_res)
np_anova_elevRange_CU1000_allAngios_res
# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2725.3, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_allAngios_res = pairwise.wilcox.test(all_angiosperms_CU1000_elevRange_coll_group_region_res$res_elevation,
                                                                    all_angiosperms_CU1000_elevRange_coll_group_region_res$region, 
                                                                    p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_allAngios_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_elevRange_coll_group_region_res$res_elevation and all_angiosperms_CU1000_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


# All angiosperms with EUDICOTS----

# Parametric Tests----
anova_elevRange_CU1000_allAngios_onlyEudicots_res = aov(res_elevation ~ region,
                                                        data =                                       all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res)
res_anova_elevRange_CU1000_allAngios_onlyEudicots_res = Anova(anova_elevRange_CU1000_allAngios_onlyEudicots_res, type = 2)
res_anova_elevRange_CU1000_allAngios_onlyEudicots_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq    Df F value    Pr(>F)    
# region     514962105     2  952.82 < 2.2e-16 ***
#   Residuals 3207089632 11868                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_elevRange_CU1000_allAngios_onlyEudicots_res = TukeyHSD(anova_elevRange_CU1000_allAngios_onlyEudicots_res)
tukey_elevRange_CU1000_allAngios_onlyEudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -455.10374 -482.12575 -428.0817     0
# Tropics-N. Temperate        97.77201   67.45564  128.0884     0
# Tropics-S. Temperate       552.87576  517.98881  587.7627     0

# Non-parametric tests----
np_anova_elevRange_CU1000_allAngios_onlyEudicots_res = kruskal.test(res_elevation ~ region,
                                                                    data = all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res)
np_anova_elevRange_CU1000_allAngios_onlyEudicots_res

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2683.5, df = 2, p-value < 2.2e-16


wrst_np_anova_elevRange_CU1000_allAngios_onlyEudicots_res = pairwise.wilcox.test(all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res$res_elevation,
                                                                                 all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res$region, 
                                                                                 p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_allAngios_onlyEudicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res$res_elevation and all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


### Test by Taxonomic Group ####

# Monocots----

# Parametric Tests----
anova_elevRange_CU1000_monocots_res = aov(res_elevation ~ region,
                                          data = all_monocots_CU1000_elevRange_coll_group_region_res)

res_anova_elevRange_CU1000_monocots_res = Anova(anova_elevRange_CU1000_monocots_res, type = 2)
res_anova_elevRange_CU1000_monocots_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region    130141680    2  211.53 < 2.2e-16 ***
#   Residuals 677070078 2201                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_monocots_res = TukeyHSD(anova_elevRange_CU1000_monocots_res)
tukey_anova_elevRange_CU1000_monocots_res

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_monocots_CU1000_elevRange_coll_group_region_res)
# 
# $region
# diff        lwr        upr     p adj
# S. Temperate-N. Temperate -526.82144 -588.96916 -464.67371 0.0000000
# Tropics-N. Temperate       -17.86906  -99.59421   63.85609 0.8651209
# Tropics-S. Temperate       508.95238  421.24369  596.66106 0.0000000


# Non-parametric tests----
np_anova_elevRange_CU1000_monocots_res = kruskal.test(res_elevation ~ region,
                                                      data = all_monocots_CU1000_elevRange_coll_group_region_res)
np_anova_elevRange_CU1000_monocots_res 

# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 603.74, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_monocots_res = pairwise.wilcox.test(all_monocots_CU1000_elevRange_coll_group_region_res$res_elevation,
                                                                   all_monocots_CU1000_elevRange_coll_group_region_res$region, 
                                                                   p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_monocots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU1000_elevRange_coll_group_region_res$res_elevation and all_monocots_CU1000_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      0.67         <2e-16      
# 
# P value adjustment method: bonferroni 


# Dicots----

# Parametric Tests----
anova_elevRange_CU1000_dicots_res = aov(res_elevation ~ region,
                                        data = all_dicots_CU1000_elevRange_coll_group_region_res)

res_anova_elevRange_CU1000_dicots_res = Anova(anova_elevRange_CU1000_dicots_res, type = 2)
res_anova_elevRange_CU1000_dicots_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     411936643    2  796.99 < 2.2e-16 ***
#   Residuals 2533432224 9803                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_dicots_res = TukeyHSD(anova_elevRange_CU1000_dicots_res)
tukey_anova_elevRange_CU1000_dicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_dicots_CU1000_elevRange_coll_group_region_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -447.1669 -476.85343 -417.4803     0
# Tropics-N. Temperate       117.3095   85.51106  149.1079     0
# Tropics-S. Temperate       564.4764  527.14339  601.8093     0

# Non-parametric tests----
np_anova_elevRange_CU1000_dicots_res = kruskal.test(res_elevation ~ region,
                                                    data = all_dicots_CU1000_elevRange_coll_group_region_res)
np_anova_elevRange_CU1000_dicots_res 
# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2161.2, df = 2, p-value < 2.2e-16

wrst_np_anova_elevRange_CU1000_dicots_res = 
  pairwise.wilcox.test(all_dicots_CU1000_elevRange_coll_group_region_res$res_elevation,
                       all_dicots_CU1000_elevRange_coll_group_region_res$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_elevRange_CU1000_dicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_elevRange_coll_group_region_res$res_elevation and all_dicots_CU1000_elevRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 

# Eudicots----

# Parametric Tests----
anova_elevRange_CU1000_eudicots_res = aov(res_elevation ~ region,
                                          data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res)

res_anova_elevRange_CU1000_eudicots_res = Anova(anova_elevRange_CU1000_eudicots_res, type = 2)
res_anova_elevRange_CU1000_eudicots_res
# Anova Table (Type II tests)
# 
# Response: res_elevation
# Sum Sq   Df F value    Pr(>F)    
# region     401074638    2  772.51 < 2.2e-16 ***
#   Residuals 2508711850 9664                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_elevRange_CU1000_eudicots_res = TukeyHSD(anova_elevRange_CU1000_eudicots_res)
tukey_anova_elevRange_CU1000_eudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_elevation ~ region, data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr       upr p adj
# S. Temperate-N. Temperate -444.5158 -474.53436 -414.4973     0
# Tropics-N. Temperate       122.2392   89.86883  154.6097     0
# Tropics-S. Temperate       566.7551  528.71591  604.7942     0


# Non-parametric tests----
np_anova_elevRange_CU1000_eudicots_res = kruskal.test(res_elevation ~ region,
                                                      data = all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res)
np_anova_elevRange_CU1000_eudicots_res 
# Kruskal-Wallis rank sum test
# 
# data:  res_elevation by region
# Kruskal-Wallis chi-squared = 2106.7, df = 2, p-value < 2.2e-16
wrst_anova_elevRange_CU1000_eudicots_res = 
  pairwise.wilcox.test(all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res$res_elevation,
                       all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res$region, 
                       p.adjust.method = "bonf")
wrst_anova_elevRange_CU1000_eudicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res$res_elevation and all_dicots_CU1000_elevRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       <2e-16      
# 
# P value adjustment method: bonferroni 


# 1.2 Differences in thermal ranges ----

# * . * . * . * . * #
# Temperature range #
# * . * . * . * . * #


### All angiosperms ####

# Parametric Tests----
anova_tempRange_precipRange_CU1000_allAngios = aov(temp_range ~ region,
                                                   data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU1000_allAngios = Anova(anova_tempRange_precipRange_CU1000_allAngios, type = 2)
res_anova_tempRange_precipRange_CU1000_allAngios
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq    Df F value    Pr(>F)    
# region    478268     2  6729.6 < 2.2e-16 ***
#   Residuals 426664 12007                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_allAngios = TukeyHSD(anova_tempRange_precipRange_CU1000_allAngios)
tukey_anova_tempRange_precipRange_CU1000_allAngios
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.258563 -12.566378 -11.9507481     0
# Tropics-N. Temperate      -13.380770 -13.723812 -13.0377277     0
# Tropics-S. Temperate       -1.122207  -1.516792  -0.7276216     0

#Non-parametric Tests----
np_anova_tempRange_precipRange_CU1000_allAngios = kruskal.test(temp_range ~ region,
                                                               data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU1000_allAngios
# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 7333.1, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_allAngios = pairwise.wilcox.test(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region$temp_range,
                                                                            all_angiosperms_CU1000_tempRange_precipRange_coll_group_region$region, 
                                                                            p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_allAngios
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_tempRange_precipRange_coll_group_region$temp_range and all_angiosperms_CU1000_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      1.2e-08     
# 
# P value adjustment method: bonferroni 

# All angiosperms with EUDICOTS----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots = aov(temp_range ~ region,
                                                                data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)

res_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots = Anova(anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots, type = 2)
res_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq    Df F value    Pr(>F)    
# region    465173     2  6528.7 < 2.2e-16 ***
#   Residuals 422802 11868                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots = TukeyHSD(anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots)
tukey_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.201722 -12.511985 -11.8914586     0
# Tropics-N. Temperate      -13.283550 -13.631639 -12.9354613     0
# Tropics-S. Temperate       -1.081828  -1.482396  -0.6812607     0


# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots = kruskal.test(temp_range ~ region,
                                                                            data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)
np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots
# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 7197.8, df = 2, p-value < 2.2e-16
wrst_np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots = pairwise.wilcox.test(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range,
                                                                                         all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$region, 
                                                                                         p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range and all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       1e-07       
# 
#P value adjustment method: bonferroni


### Test by Taxonomic Group ####

# Monocots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_monocots = aov(temp_range ~ region,
                                                  data = all_monocots_CU1000_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU1000_monocots = Anova(anova_tempRange_precipRange_CU1000_monocots, type = 2)
res_anova_tempRange_precipRange_CU1000_monocots
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region     96726    2  1346.6 < 2.2e-16 ***
#   Residuals  79048 2201                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_monocots = TukeyHSD(anova_tempRange_precipRange_CU1000_monocots)
tukey_anova_tempRange_precipRange_CU1000_monocots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_monocots_CU1000_tempRange_precipRange_coll_group_region)
# 
# $region
# diff        lwr         upr     p adj
# S. Temperate-N. Temperate -12.757509 -13.429022 -12.0859967 0.0000000
# Tropics-N. Temperate      -14.320661 -15.203710 -13.4376124 0.0000000
# Tropics-S. Temperate       -1.563152  -2.510853  -0.6154504 0.0003318

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_monocots = kruskal.test(temp_range ~ region,
                                                              data = all_monocots_CU1000_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU1000_monocots 

# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 1408, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_monocots  = 
  pairwise.wilcox.test(all_monocots_CU1000_tempRange_precipRange_coll_group_region$temp_range,
                       all_monocots_CU1000_tempRange_precipRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_monocots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU1000_tempRange_precipRange_coll_group_region$temp_range and all_monocots_CU1000_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      1.2e-05     
# 
# P value adjustment method: bonferroni 

# Dicots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_dicots = aov(temp_range ~ region,
                                                data = all_dicots_CU1000_tempRange_precipRange_coll_group_region)

res_anova_tempRange_precipRange_CU1000_dicots = Anova(anova_tempRange_precipRange_CU1000_dicots, type = 2)
res_anova_tempRange_precipRange_CU1000_dicots
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    382983    2  5428.1 < 2.2e-16 ***
#   Residuals 345826 9803                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_dicots = TukeyHSD(anova_tempRange_precipRange_CU1000_dicots)
tukey_anova_tempRange_precipRange_CU1000_dicots
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_dicots_CU1000_tempRange_precipRange_coll_group_region)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.2079134 -12.554757 -11.8610696 0e+00
# Tropics-N. Temperate      -13.1814047 -13.552923 -12.8098867 0e+00
# Tropics-S. Temperate       -0.9734913  -1.409672  -0.5373103 5e-07

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_dicots = kruskal.test(temp_range ~ region,
                                                            data = all_dicots_CU1000_tempRange_precipRange_coll_group_region)
np_anova_tempRange_precipRange_CU1000_dicots 
# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 5932.7, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_dicots = 
  pairwise.wilcox.test(all_dicots_CU1000_tempRange_precipRange_coll_group_region$temp_range,
                       all_dicots_CU1000_tempRange_precipRange_coll_group_region$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_dicots

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_tempRange_precipRange_coll_group_region$temp_range and all_dicots_CU1000_tempRange_precipRange_coll_group_region$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      0.00025     
# 
# P value adjustment method: bonferroni 


# Eudicots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_eudicots = aov(temp_range ~ region,
                                                  data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)

res_anova_tempRange_precipRange_CU1000_eudicots = Anova(anova_tempRange_precipRange_CU1000_eudicots, type = 2)
res_anova_tempRange_precipRange_CU1000_eudicots
# Anova Table (Type II tests)
# 
# Response: temp_range
# Sum Sq   Df F value    Pr(>F)    
# region    370005    2  5227.3 < 2.2e-16 ***
#   Residuals 342024 9664                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_eudicots = TukeyHSD(anova_tempRange_precipRange_CU1000_eudicots)
tukey_anova_tempRange_precipRange_CU1000_eudicots 
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temp_range ~ region, data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)
# 
# $region
# diff        lwr         upr   p adj
# S. Temperate-N. Temperate -12.1364107 -12.486914 -11.7859071 0.0e+00
# Tropics-N. Temperate      -13.0669543 -13.444919 -12.6889895 0.0e+00
# Tropics-S. Temperate       -0.9305436  -1.374698  -0.4863892 2.8e-06

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_eudicots = kruskal.test(temp_range ~ region,
                                                              data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots)
np_anova_tempRange_precipRange_CU1000_eudicots 
# Kruskal-Wallis rank sum test
# 
# data:  temp_range by region
# Kruskal-Wallis chi-squared = 5795, df = 2, p-value < 2.2e-16

wrst_anova_tempRange_precipRange_CU1000_eudicots = 
  pairwise.wilcox.test(all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range,
                       all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$region, 
                       p.adjust.method = "bonf")
wrst_anova_tempRange_precipRange_CU1000_eudicots
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$temp_range and all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       0.001       
# 
# P value adjustment method: bonferroni  


### Analysis with residuals ####

## Perform all the tests again, but this time use residuals of linear regression (elevation_range ~ n_collections) as the response variable.

# First estimate residuals, test for normality./parametric test, then run tests.

# Estimate residuals----

###  All Angiosperms ####

# Elevation range----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region %>% 
  ungroup() %>%
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

#Plot to check
# all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res %>%
#   ggplot(aes(x = region,
#              y = res_temp)) +
#   geom_boxplot() +
#   theme_classic()

# Testing for parametric analysis ----

all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
# group     2  411.69 < 2.2e-16 ***
#       12007                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Only with Eudicots ----


# UNCOMMENT BELOW FOR RESIDUALS PER REGION!
# all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()


# These are the "overall residuals", not the residuals per region.
all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

# Testing for parametric analysis ----


all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# Levene's Test for Homogeneity of Variance (center = median)
#          Df F value    Pr(>F)    
# group     2  411.92 < 2.2e-16 ***
#       11868                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###  By Taxonomic groups ####

# Monocots----

#UNCOMMENT FR RESIDUASL PER REGION
# all_monocots_CU1000_tempRange_precipRange_coll_group_region_res = all_monocots_CU1000_tempRange_precipRange_coll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the "overall" residuals
all_monocots_CU1000_tempRange_precipRange_coll_group_region_res = all_monocots_CU1000_tempRange_precipRange_coll_group_region %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))

# Testing for parametric analysis ----

all_monocots_CU1000_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_monocots_CU1000_tempRange_precipRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  121.37 < 2.2e-16 ***
#       2201                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Dicots----

# UNCOMMNET BELOW FOR RESIDUALS PER REGION
# all_dicots_CU1000_tempRange_precipRange_coll_group_region_res = all_dicots_CU1000_tempRange_precipRangecoll_group_region %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_elevation = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()


# These are the overall residuals
all_dicots_CU1000_tempRange_precipRange_coll_group_region_res = all_dicots_CU1000_tempRange_precipRange_coll_group_region %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))


# Testing for parametric analysis ----

all_dicots_CU1000_tempRange_precipRange_coll_group_region_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_dicots_CU1000_tempRange_precipRange_coll_group_region_res)
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2  298.13 < 2.2e-16 ***
#       9803                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Eudicots----

# UNCOMMENT BELOW FOR RESIDUALS PER REGION
# all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
#   ungroup() %>%
#   group_by(region) %>%
#   mutate(res_temp = residuals(lm(temp_range ~ n_collections))) %>% 
#   ungroup()

# These are the overall residuals
all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots %>% 
  mutate(res_temp = residuals(lm(temp_range ~ n_collections)))


# Testing for parametric analysis ----

all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res %>% 
  ggplot(aes(sample = res_temp)) +
  stat_qq() +
  stat_qq_line() + 
  theme_classic()

# test for Homoscedasticity (homogeneity of variance)
car::leveneTest(res_temp ~ region, all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
# group    2     298 < 2.2e-16 ***
#       9664                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA and KW with residuals----

# Parametric Tests----
# All angiosperms----
anova_tempRange_precipRange_CU1000_allAngios_res = aov(res_temp ~ region,
                                                       data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU1000_allAngios_res = Anova(anova_tempRange_precipRange_CU1000_allAngios_res, type = 2)
res_anova_tempRange_precipRange_CU1000_allAngios_res
# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq    Df F value    Pr(>F)    
# region    477417     2  6704.9 < 2.2e-16 ***
#   Residuals 427476 12007                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_allAngios_res = TukeyHSD(anova_tempRange_precipRange_CU1000_allAngios_res)
tukey_anova_tempRange_precipRange_CU1000_allAngios_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.230853 -12.538961 -11.9227457     0
# Tropics-N. Temperate      -13.388216 -13.731584 -13.0448472     0
# Tropics-S. Temperate       -1.157362  -1.552323  -0.7624018     0

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_allAngios_res = kruskal.test(res_temp ~ region,
                                                                   data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU1000_allAngios_res
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 7316.7, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_allAngios_res = pairwise.wilcox.test(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res$res_temp,
                                                                                all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res$region, 
                                                                                p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_allAngios_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res$res_temp and all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      9.5e-10     
# 
# P value adjustment method: bonferroni  

# All angiosperms with EUDICOTS----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res = aov(res_temp ~ region,
                                                                    data =                                       all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
res_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res = Anova(anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res, type = 2)
res_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res
# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq    Df F value    Pr(>F)    
# region    464474     2  6508.5 < 2.2e-16 ***
#   Residuals 423475 11868                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res = TukeyHSD(anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res)
tukey_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.178756 -12.489266 -11.8682453     0
# Tropics-N. Temperate      -13.289808 -13.638174 -12.9414423     0
# Tropics-S. Temperate       -1.111053  -1.511939  -0.7101661     0

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res = kruskal.test(res_temp ~ region,
                                                                                data = all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 7183.1, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res = pairwise.wilcox.test(all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp,
                                                                                             all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$region, 
                                                                                             p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_allAngios_onlyEudicots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp and all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      1.2e-08     
# 
# P value adjustment method: bonferroni  


### Test by Taxonomic Group ####

# Monocots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_monocots_res = aov(res_temp ~ region,
                                                      data = all_monocots_CU1000_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU1000_monocots_res = Anova(anova_tempRange_precipRange_CU1000_monocots_res, type = 2)
res_anova_tempRange_precipRange_CU1000_monocots_res
# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region     97457    2  1373.5 < 2.2e-16 ***
#   Residuals  78089 2201                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_monocots_res = TukeyHSD(anova_tempRange_precipRange_CU1000_monocots_res)
tukey_anova_tempRange_precipRange_CU1000_monocots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_monocots_CU1000_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff        lwr         upr     p adj
# S. Temperate-N. Temperate -12.878665 -13.546091 -12.2112386 0.0000000
# Tropics-N. Temperate      -14.245694 -15.123370 -13.3680191 0.0000000
# Tropics-S. Temperate       -1.367029  -2.308964  -0.4250948 0.0019555

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_monocots_res = kruskal.test(res_temp ~ region,
                                                                  data = all_monocots_CU1000_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU1000_monocots_res 
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 1420.1, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_monocots_res = pairwise.wilcox.test(all_monocots_CU1000_tempRange_precipRange_coll_group_region_res$res_temp,
                                                                               all_monocots_CU1000_tempRange_precipRange_coll_group_region_res$region, 
                                                                               p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_monocots_res
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_monocots_CU1000_tempRange_precipRange_coll_group_region_res$res_temp and all_monocots_CU1000_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate <2e-16       -           
#   Tropics      <2e-16       2e-04       
# 
# P value adjustment method: bonferroni 


# Dicots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_dicots_res = aov(res_temp ~ region,
                                                    data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_res)

res_anova_tempRange_precipRange_CU1000_dicots_res = Anova(anova_tempRange_precipRange_CU1000_dicots_res, type = 2)
res_anova_tempRange_precipRange_CU1000_dicots_res
# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    381405    2  5383.6 < 2.2e-16 ***
#   Residuals 347249 9803                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_dicots_res = TukeyHSD(anova_tempRange_precipRange_CU1000_dicots_res)
tukey_anova_tempRange_precipRange_CU1000_dicots_res
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_res)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.143251 -12.490808 -11.7956941 0e+00
# Tropics-N. Temperate      -13.195706 -13.567987 -12.8234241 0e+00
# Tropics-S. Temperate       -1.052455  -1.489532  -0.6153773 1e-07

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_dicots_res = kruskal.test(res_temp ~ region,
                                                                data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_res)
np_anova_tempRange_precipRange_CU1000_dicots_res 
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 5905.8, df = 2, p-value < 2.2e-16

wrst_np_anova_tempRange_precipRange_CU1000_dicots_res = 
  pairwise.wilcox.test(all_dicots_CU1000_tempRange_precipRange_coll_group_region_res$res_temp,
                       all_dicots_CU1000_tempRange_precipRange_coll_group_region_res$region, 
                       p.adjust.method = "bonf")
wrst_np_anova_tempRange_precipRange_CU1000_dicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_tempRange_precipRange_coll_group_region_res$res_temp and all_dicots_CU1000_tempRange_precipRange_coll_group_region_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      1.4e-05     
# 
# P value adjustment method: bonferroni 

# Eudicots----

# Parametric Tests----
anova_tempRange_precipRange_CU1000_eudicots_res = aov(res_temp ~ region,
                                                      data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)

res_anova_tempRange_precipRange_CU1000_eudicots_res = Anova(anova_tempRange_precipRange_CU1000_eudicots_res, type = 2)
res_anova_tempRange_precipRange_CU1000_eudicots_res
# Anova Table (Type II tests)
# 
# Response: res_temp
# Sum Sq   Df F value    Pr(>F)    
# region    368556    2  5186.8 < 2.2e-16 ***
#   Residuals 343343 9664                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_anova_tempRange_precipRange_CU1000_eudicots_res = TukeyHSD(anova_tempRange_precipRange_CU1000_eudicots_res)
tukey_anova_tempRange_precipRange_CU1000_eudicots_res

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = res_temp ~ region, data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
# 
# $region
# diff        lwr         upr p adj
# S. Temperate-N. Temperate -12.076298 -12.427476 -11.7251190 0e+00
# Tropics-N. Temperate      -13.080434 -13.459127 -12.7017417 0e+00
# Tropics-S. Temperate       -1.004137  -1.449147  -0.5591269 4e-07

# Non-parametric tests----
np_anova_tempRange_precipRange_CU1000_eudicots_res = kruskal.test(res_temp ~ region,
                                                                  data = all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res)
np_anova_tempRange_precipRange_CU1000_eudicots_res 
# Kruskal-Wallis rank sum test
# 
# data:  res_temp by region
# Kruskal-Wallis chi-squared = 5769.4, df = 2, p-value < 2.2e-16

wrst_anova_tempRange_precipRange_CU1000_eudicots_res = 
  pairwise.wilcox.test(all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp,
                       all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$region, 
                       p.adjust.method = "bonf")
wrst_anova_tempRange_precipRange_CU1000_eudicots_res

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$res_temp and all_dicots_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots_res$region 
# 
# N. Temperate S. Temperate
# S. Temperate < 2e-16      -           
#   Tropics      < 2e-16      8.7e-05     
# 
# P value adjustment method: bonferroni 


#-----------------------------------------------------------------------#
#-----------------------END OF ANALYSES PREDICTION 1--------------------#
#-----------------------------------------------------------------------#



#### PREDICTION 2: Differences in Overlap #### 

#<><><><><><><><><><><><>#
#<><><><><><><><><><><><>#
#                        #
# DIFFERENCES IN OVERLAP #
#                        #
#<><><><><><><><><><><><>#
#<><><><><><><><><><><><>#


# CU 500 ----

# Elevation Overlap  ====

# All angiosperms ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
#total_n_collections = c(species_pairs_elevation_range_CU500_wide_ElevOverlap$sp1_n_collections + species_pairs_elevation_range_CU500_wide_ElevOverlap$sp2_n_collections)


species_pairs_elevation_range_CU500_wide_ElevOverlap %>% 
  mutate(total_n_collections = sp1_n_collections + sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU500_int_2cov = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap)

Anova(ancovaPar_elevOverlap_CU500_int_2cov, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df  F value  Pr(>F)    
# (Intercept)         67.583   1 984.2857 < 2e-16 ***
#   total_n_collections  0.008   1   0.1174 0.73213    
# pair_age             0.040   1   0.5762 0.44844    
# sp1_region           0.177   2   1.2877 0.27751    
# pair_age:sp1_region  0.350   2   2.5504 0.07984 .  
# Residuals           19.500 284                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU500_Noint_2cov = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap)
Anova(ancovaPar_elevOverlap_CU500_Noint_2cov, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value   Pr(>F)   
# total_n_collections  0.0053   1  0.0761 0.782807   
# pair_age             0.1208   1  1.7400 0.188196   
# sp1_region           0.8756   2  6.3075 0.002086 **
#   Residuals           19.8503 286                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Region alone is significant 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_elevOverlap_CU500_Noint_2cov, "sp1_region", ordered  =TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = elevation_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_elevation_range_CU500_wide_ElevOverlap)
# 
# $sp1_region
# diff          lwr       upr     p adj
# S. Temperate-N. Temperate 0.10067877  0.023111856 0.1782457 0.0068711
# Tropics-N. Temperate      0.11541159 -0.003599548 0.2344227 0.0595673
# Tropics-S. Temperate      0.01473283 -0.105744625 0.1352103 0.9552830


# Non parametric ANCOVA ----
# NOTE: NOT ABLE TO RUN THIS YET



# All angiosperms with Eudicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots$sp1_n_collections + species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU500_int_2cov_onlyEudicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots)

Anova(ancovaPar_elevOverlap_CU500_int_2cov_onlyEudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df  F value  Pr(>F)    
# (Intercept)         67.469   1 979.3215 < 2e-16 ***
#   total_n_collections  0.008   1   0.1146 0.73519    
# pair_age             0.041   1   0.5956 0.44090    
# sp1_region           0.180   2   1.3044 0.27296    
# pair_age:sp1_region  0.336   2   2.4378 0.08919 .  
# Residuals           19.497 283                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots)
Anova(ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value   Pr(>F)   
# total_n_collections  0.0051   1  0.0735 0.786491   
# pair_age             0.1383   1  1.9875 0.159695   
# sp1_region           0.8786   2  6.3132 0.002076 **
#   Residuals           19.8327 285                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Region alone is significant 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots, "sp1_region", ordered  =TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = elevation_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots)
# 
# $sp1_region
# diff          lwr       upr     p adj
# S. Temperate-N. Temperate 0.10113373  0.023292747 0.1789747 0.0068070
# Tropics-N. Temperate      0.11591898 -0.003250212 0.2350882 0.0585640
# Tropics-S. Temperate      0.01478524 -0.105962416 0.1355329 0.9551667



# Non parametric ANCOVA ----




# Test by taxonomic groups ----

# Monocots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots$sp1_n_collections + species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU500_int_2cov_monocots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots)

Anova(ancovaPar_elevOverlap_CU500_int_2cov_monocots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq Df  F value    Pr(>F)    
# (Intercept)         13.0023  1 135.2882 5.831e-16 ***
#   total_n_collections  0.1880  1   1.9564    0.1680    
# pair_age             0.0830  1   0.8631    0.3572    
# sp1_region           0.0495  2   0.2576    0.7739    
# pair_age:sp1_region  0.0498  2   0.2591    0.7727    
# Residuals            4.9015 51                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU500_Noint_2cov_monocots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots)
Anova(ancovaPar_elevOverlap_CU500_Noint_2cov_monocots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq Df F value Pr(>F)
# total_n_collections 0.2032  1  2.1753 0.1462
# pair_age            0.2151  1  2.3027 0.1351
# sp1_region          0.2636  2  1.4107 0.2530
# Residuals           4.9513 53   

# No covariate is significant


# Non parametric ANCOVA ----




# Dicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots$sp1_n_collections + species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU500_int_2cov_dicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots)

Anova(ancovaPar_elevOverlap_CU500_int_2cov_dicots, type = 3 )

# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df  F value Pr(>F)    
# (Intercept)         52.541   1 832.5431 <2e-16 ***
#   total_n_collections  0.037   1   0.5895 0.4434    
# pair_age             0.008   1   0.1316 0.7172    
# sp1_region           0.195   2   1.5458 0.2154    
# pair_age:sp1_region  0.234   2   1.8544 0.1589    
# Residuals           14.263 226                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU500_Noint_2cov_dicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots)
Anova(ancovaPar_elevOverlap_CU500_Noint_2cov_dicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value   Pr(>F)   
# total_n_collections  0.0328   1  0.5161 0.473249   
# pair_age             0.0624   1  0.9817 0.322840   
# sp1_region           0.7221   2  5.6781 0.003922 **
#   Residuals           14.4967 228                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_elevOverlap_CU500_Noint_2cov_dicots, "sp1_region", ordered  =TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = elevation_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots)
# 
# $sp1_region
# diff         lwr       upr     p adj
# S. Temperate-N. Temperate 0.10293235  0.02023953 0.1856252 0.0101946
# Tropics-N. Temperate      0.11698946 -0.01344317 0.2474221 0.0888456
# Tropics-S. Temperate      0.01405712 -0.11650208 0.1446163 0.9650648


# Non parametric ANCOVA ----


# Eudicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots$sp1_n_collections + species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU500_int_2cov_onlyEudicots_eudicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots)

Anova(ancovaPar_elevOverlap_CU500_int_2cov_onlyEudicots_eudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df  F value Pr(>F)    
# (Intercept)         52.467   1 827.8581 <2e-16 ***
#   total_n_collections  0.037   1   0.5819 0.4464    
# pair_age             0.009   1   0.1412 0.7074    
# sp1_region           0.198   2   1.5597 0.2125    
# pair_age:sp1_region  0.227   2   1.7874 0.1698    
# Residuals           14.260 225                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots_eudicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots)
Anova(ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots_eudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value  Pr(>F)   
# total_n_collections  0.0325   1  0.5090 0.47632   
# pair_age             0.0726   1  1.1380 0.28721   
# sp1_region           0.7245   2  5.6764 0.00393 **
#   Residuals           14.4865 227                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_elevOverlap_CU500_Noint_2cov_onlyEudicots_eudicots, "sp1_region", ordered  =TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = elevation_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_elevation_range_CU500_wide_ElevOverlap_onlyEudicots_eudicots)
# 
# $sp1_region
# diff         lwr       upr     p adj
# S. Temperate-N. Temperate 0.10346073  0.02040910 0.1865124 0.0101247
# Tropics-N. Temperate      0.11744278 -0.01323437 0.2481199 0.0880037
# Tropics-S. Temperate      0.01398205 -0.11695109 0.1449152 0.9656221


# Non parametric ANCOVA ----



# Temperature Overlap ----

# Parametric ANCOVA ----

# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU500_wide_TempOverlap$sp1_n_collections + species_pairs_bioclim_CU500_wide_TempOverlap$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU500_int_2cov = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap)

Anova(ancovaPar_TempOverlap_CU500_int_2cov, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         9014.1   1 2045.7490 < 2.2e-16 ***
#   total_n_collections    1.5   1    0.3503 0.5544209    
# pair_age              10.9   1    2.4812 0.1163254    
# sp1_region            68.5   2    7.7758 0.0005155 ***
#   pair_age:sp1_region    6.3   2    0.7096 0.4927085    
# Residuals           1251.4 284                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_TempOverlap_CU500_Noint_2cov = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap)
Anova(ancovaPar_TempOverlap_CU500_Noint_2cov, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections    1.78   1  0.4052    0.5249    
# pair_age              24.98   1  5.6809    0.0178 *  
#   sp1_region           143.73   2 16.3430 1.903e-07 ***
#   Residuals           1257.63 286                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Both pair age and region are significant!

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU500_Noint_2cov, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU500_wide_TempOverlap)
# 
# $sp1_region
# diff        lwr      upr     p adj
# S. Temperate-N. Temperate 1.2807911  0.6633872 1.898195 0.0000051
# Tropics-N. Temperate      1.5094810  0.5621962 2.456766 0.0006154
# Tropics-S. Temperate      0.2286898 -0.7302662 1.187646 0.8404493

# Non parametric ANCOVA ----



# All angiosperms with Eudicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots$sp1_n_collections + species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU500_int_2cov_onlyEudicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots)

Anova(ancovaPar_TempOverlap_CU500_int_2cov_onlyEudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         8976.1   1 2035.2690 < 2.2e-16 ***
#   total_n_collections    1.5   1    0.3338 0.5638691    
# pair_age              10.0   1    2.2688 0.1331143    
# sp1_region            65.8   2    7.4549 0.0006995 ***
#   pair_age:sp1_region    7.6   2    0.8641 0.4225264    
# Residuals           1248.1 283                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots)
Anova(ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value   Pr(>F)    
# total_n_collections    1.75   1  0.3971  0.52910    
# pair_age              17.79   1  4.0378  0.04543 *  
#   sp1_region           143.27   2 16.2579 2.06e-07 ***
#   Residuals           1255.73 285                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots, "sp1_region", ordered  = FALSE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots)
# 
# $sp1_region
# diff        lwr      upr     p adj
# S. Temperate-N. Temperate 1.2845715  0.6651788 1.903964 0.0000051
# Tropics-N. Temperate      1.5030071  0.5547596 2.451255 0.0006635
# Tropics-S. Temperate      0.2184356 -0.7423720 1.179243 0.8538552


# Non parametric ANCOVA ----




# Test by taxonomic groups ----

# Monocots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU500_wide_TempOverlap_monocots$sp1_n_collections + species_pairs_bioclim_CU500_wide_TempOverlap_monocots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU500_int_2cov_monocots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_monocots)

Anova(ancovaPar_TempOverlap_CU500_int_2cov_monocots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq Df  F value  Pr(>F)    
# (Intercept)         1495.72  1 320.7565 < 2e-16 ***
#   total_n_collections    0.90  1   0.1932 0.66211    
# pair_age              12.85  1   2.7568 0.10298    
# sp1_region            34.21  2   3.6679 0.03249 *  
#   pair_age:sp1_region   40.16  2   4.3057 0.01871 *  
#   Residuals            237.82 51                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU500_int_2cov_monocots, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend    SE df lower.CL upper.CL
# N. Temperate         -0.224 0.169 51  -0.5630   0.1155
# S. Temperate          0.338 0.214 51  -0.0913   0.7667
# Tropics              -0.986 0.449 51  -1.8887  -0.0843
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate    SE df t.ratio p.value
# N. Temperate - S. Temperate   -0.561 0.272 51  -2.062  0.1081
# N. Temperate - Tropics         0.763 0.480 51   1.589  0.2598
# S. Temperate - Tropics         1.324 0.498 51   2.661  0.0276
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# We see ths slopes and that the slope ST vs Tr is the only significantly different. I believe this islargely driven because of sample size, ie. the tropics have very low sample size. NT and Tr are not significantly different because both slopes are negative.

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU500_int_2cov_monocots, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))
# The cov.reduce = range argument is passed to ref_grid(); it is needed because by default, each covariate is reduced to only one value, in our case the average of pair_ages. The current function cannot accommodate all values. See https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html Instead we use the range for visualization

# This result means that in the Tropics and N.Temp as species divergence increases, the termal overlap decreases, whereas the opposite is true in the S. Temp. In other words, in the Tropics and N. Temp recently diverged species show more overlap and with time they "evolve" dis-similarities in thermal overlap.

# Because age was significant, run a model without interactions
ancovaPar_TempOverlap_CU500_Noint_2cov_monocots = aov(temperature_overlap ~ total_n_collections + pair_age + sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_monocots)

Anova(ancovaPar_TempOverlap_CU500_Noint_2cov_monocots, type = 2 )
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq Df F value Pr(>F)
# total_n_collections   0.257  1  0.0489 0.8258
# pair_age              2.144  1  0.4087 0.5254
# sp1_region           13.178  2  1.2562 0.2931
# Residuals           277.973 53 

# No significant factors

# Non parametric ANCOVA ----





# Dicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU500_wide_TempOverlap_dicots$sp1_n_collections + species_pairs_bioclim_CU500_wide_TempOverlap_dicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU500_int_2cov_dicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_dicots)

Anova(ancovaPar_TempOverlap_CU500_int_2cov_dicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         7112.4   1 1671.3494 < 2.2e-16 ***
#   total_n_collections    2.5   1    0.5811 0.4466899    
# pair_age               2.5   1    0.5862 0.4447079    
# sp1_region            76.2   2    8.9497 0.0001817 ***
#   pair_age:sp1_region    5.0   2    0.5822 0.5594887    
# Residuals            961.7 226                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Interaction is non-significant. Run a model without interaction
ancovaPar_TempOverlap_CU500_Noint_2cov_dicots = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_dicots)
Anova(ancovaPar_TempOverlap_CU500_Noint_2cov_dicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections   2.56   1  0.6039   0.43790    
# pair_age             25.93   1  6.1154   0.01413 *  
#   sp1_region          139.87   2 16.4941 2.041e-07 ***
#   Residuals           966.69 228                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Both age and region are significant. 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU500_Noint_2cov_dicots, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU500_wide_TempOverlap_dicots)
# 
# $sp1_region
# diff        lwr      upr     p adj
# S. Temperate-N. Temperate 1.45274294  0.7774728 2.128013 0.0000024
# Tropics-N. Temperate      1.54393799  0.4788243 2.609052 0.0021371
# Tropics-S. Temperate      0.09119505 -0.9749522 1.157342 0.9778040


# Non parametric ANCOVA ----



# Eudicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots$sp1_n_collections + species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU500_int_2cov_onlyEudicots_eudicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots)

Anova(ancovaPar_TempOverlap_CU500_int_2cov_onlyEudicots_eudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         7085.8   1 1661.9489 < 2.2e-16 ***
#   total_n_collections    2.4   1    0.5619 0.4542872    
# pair_age               2.1   1    0.5034 0.4787383    
# sp1_region            73.3   2    8.5995 0.0002518 ***
#   pair_age:sp1_region    5.4   2    0.6274 0.5348945    
# Residuals            959.3 225                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction not significant. Run mdoel without interactions
ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots_eudicots = aov(temperature_overlap ~ total_n_collections + pair_age + sp1_region, species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots)

Anova(ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots_eudicots, type = 2 )

# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections   2.52   1  0.5928   0.44215    
# pair_age             18.18   1  4.2778   0.03975 *  
#   sp1_region          139.35   2 16.3958 2.233e-07 ***
#   Residuals           964.65 227                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Both age and region are significant. 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU500_Noint_2cov_onlyEudicots_eudicots, "sp1_region", ordered  = TRUE)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU500_wide_TempOverlap_onlyEudicots_eudicots)
# 
# $sp1_region
# diff        lwr      upr     p adj
# S. Temperate-N. Temperate 1.4573104  0.7795870 2.135034 0.0000024
# Tropics-N. Temperate      1.5357301  0.4693699 2.602090 0.0023068
# Tropics-S. Temperate      0.0784197 -0.9900294 1.146869 0.9836078



### CU 1000 ####

# Elevation Overlap ----

# All angiosperms ====

# Parametric ANCOVA ----

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU1000_wide_ElevOverlap$sp1_n_collections + species_pairs_elevation_range_CU1000_wide_ElevOverlap$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU1000_int_2cov = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap)

Anova(ancovaPar_elevOverlap_CU1000_int_2cov, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df   F value Pr(>F)    
# (Intercept)         91.401   1 1587.0374 <2e-16 ***
#   total_n_collections  0.104   1    1.8035 0.1799    
# pair_age             0.077   1    1.3385 0.2478    
# sp1_region           0.142   2    1.2306 0.2930    
# pair_age:sp1_region  0.148   2    1.2847 0.2776    
# Residuals           29.602 514                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU1000_Noint_2cov = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap)
Anova(ancovaPar_elevOverlap_CU1000_Noint_2cov, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value Pr(>F)
# total_n_collections  0.1162   1  2.0153 0.1563
# pair_age             0.0043   1  0.0739 0.7858
# sp1_region           0.0910   2  0.7889 0.4549
# Residuals           29.7503 516              

# No factor is significant. There are NO significant differences in elevation overlap


# Non parametric ANCOVA ----
# NOTE: NOT ABLE TO RUN THIS YET





# All angiosperms with Eudicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots$sp1_n_collections + species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU1000_int_2cov_onlyEudicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots)

Anova(ancovaPar_elevOverlap_CU1000_int_2cov_onlyEudicots, type = 3 )

# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df   F value Pr(>F)    
# (Intercept)         91.209   1 1578.1464 <2e-16 ***
#   total_n_collections  0.103   1    1.7738 0.1835    
# pair_age             0.072   1    1.2384 0.2663    
# sp1_region           0.149   2    1.2849 0.2776    
# pair_age:sp1_region  0.157   2    1.3586 0.2579    
# Residuals           29.591 512                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU1000_Noint_2cov_onlyEudicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots)
Anova(ancovaPar_elevOverlap_CU1000_Noint_2cov_onlyEudicots, type = 2)

# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value Pr(>F)
# total_n_collections  0.1162   1  2.0075 0.1571
# pair_age             0.0057   1  0.0990 0.7531
# sp1_region           0.0907   2  0.7838 0.4572
# Residuals           29.7481 514  

# No factor is significant. NO significant differences in overlap


# Non parametric ANCOVA ----





# Test by taxonomic groups ----

# Monocots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots$sp1_n_collections + species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU1000_int_2cov_monocots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots)

Anova(ancovaPar_elevOverlap_CU1000_int_2cov_monocots, type = 3 )

# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq Df  F value Pr(>F)    
# (Intercept)         18.5917  1 485.8070 <2e-16 ***
#   total_n_collections  0.0020  1   0.0534 0.8178    
# pair_age             0.0008  1   0.0197 0.8886    
# sp1_region           0.0468  2   0.6118 0.5445    
# pair_age:sp1_region  0.0373  2   0.4867 0.6162    
# Residuals            3.5973 94                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU1000_Noint_2cov_monocots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots)
Anova(ancovaPar_elevOverlap_CU1000_Noint_2cov_monocots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq Df F value Pr(>F)
# total_n_collections 0.0043  1  0.1132 0.7372
# pair_age            0.0017  1  0.0446 0.8332
# sp1_region          0.0868  2  1.1464 0.3221
# Residuals           3.6346 96  

# No factor is significant. NO differences in overlap


# Non parametric ANCOVA ----




# Dicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots$sp1_n_collections + species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU1000_int_2cov_dicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots)

Anova(ancovaPar_elevOverlap_CU1000_int_2cov_dicots, type = 3 )

# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df   F value Pr(>F)    
# (Intercept)         69.077   1 1122.0316 <2e-16 ***
#   total_n_collections  0.082   1    1.3356 0.2485    
# pair_age             0.063   1    1.0242 0.3121    
# sp1_region           0.236   2    1.9206 0.1478    
# pair_age:sp1_region  0.123   2    0.9994 0.3690    
# Residuals           25.426 413                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is not significant. Run model without interactions

ancovaPar_elevOverlap_CU1000_Noint_2cov_dicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots)
Anova(ancovaPar_elevOverlap_CU1000_Noint_2cov_dicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value Pr(>F)
# total_n_collections  0.0900   1  1.4614 0.2274
# pair_age             0.0122   1  0.1982 0.6564
# sp1_region           0.2481   2  2.0147 0.1347
# Residuals           25.5492 415 

# No factor is significant. No significant differences

# Non parametric ANCOVA ----





# Eudicots ----

# Parametric ANCOVA ----
# Run ANCOVA with interactions and 2 covariates, divergence time and collections per regions

# Get total number of samples
total_n_collections = c(species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots$sp1_n_collections + species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_elevOverlap_CU1000_int_2cov_onlyEudicots_eudicots = aov(elevation_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots)

Anova(ancovaPar_elevOverlap_CU1000_int_2cov_onlyEudicots_eudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: elevation_overlap
# Sum Sq  Df   F value Pr(>F)    
# (Intercept)         68.959   1 1115.1491 <2e-16 ***
#   total_n_collections  0.081   1    1.3117 0.2527    
# pair_age             0.059   1    0.9534 0.3294    
# sp1_region           0.245   2    1.9773 0.1398    
# pair_age:sp1_region  0.130   2    1.0523 0.3501    
# Residuals           25.416 411                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is non-significant. Run a model without interaction
ancovaPar_elevOverlap_CU1000_Noint_2cov_onlyEudicots_eudicots = aov(elevation_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_elevation_range_CU1000_wide_ElevOverlap_onlyEudicots_eudicots)
Anova(ancovaPar_elevOverlap_CU1000_Noint_2cov_onlyEudicots_eudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: elevation_overlap
# Sum Sq  Df F value Pr(>F)
# total_n_collections  0.0899   1  1.4537 0.2286
# pair_age             0.0149   1  0.2402 0.6244
# sp1_region           0.2483   2  2.0070 0.1357
# Residuals           25.5457 413 

# No factor is significant. No significant differences in overlap


# Non parametric ANCOVA ----




# Temperature Overlap ----

# All angisoperms ----

# Parametric ANCOVA ----

# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU1000_wide_TempOverlap$sp1_n_collections + species_pairs_bioclim_CU1000_wide_TempOverlap$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particualry when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU1000_int_2cov = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap)

Anova(ancovaPar_TempOverlap_CU1000_int_2cov, type = 3 )

# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         12201.7   1 2605.5026 < 2.2e-16 ***
#   total_n_collections     4.6   1    0.9908  0.320009    
# pair_age                0.8   1    0.1742  0.676574    
# sp1_region             45.5   2    4.8567  0.008136 ** 
#   pair_age:sp1_region    44.1   2    4.7124  0.009375 ** 
#   Residuals            2407.1 514                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU1000_int_2cov, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend     SE  df lower.CL upper.CL
# N. Temperate        -0.1615 0.0389 514  -0.2378  -0.0851
# S. Temperate        -0.0365 0.0304 514  -0.0962   0.0232
# Tropics              0.1410 0.1272 514  -0.1088   0.3909
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate     SE  df t.ratio p.value
# N. Temperate - S. Temperate   -0.125 0.0494 514  -2.530  0.0314
# N. Temperate - Tropics        -0.302 0.1330 514  -2.275  0.0603
# S. Temperate - Tropics        -0.178 0.1308 514  -1.358  0.3640
# 
# P value adjustment: tukey method for comparing a family of 3 estimates  

# NT vs ST is significant and NT vs Tr is almost significant (0.06).

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU1000_int_2cov, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))

# Because region is significant, let's run a model without interactions considering only regions as main effect
ancovaPar_TempOverlap_CU1000_Noint_2cov = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap)
Anova(ancovaPar_TempOverlap_CU1000_Noint_2cov, type = 2)
# Anova Table (Type II tests)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections    6.47   1  1.3614  0.243843    
# pair_age              49.35   1 10.3885  0.001348 ** 
#   sp1_region           127.33   2 13.4016 2.118e-06 ***
#   Residuals           2451.23 516                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Tukey with region as main effect
TukeyHSD(ancovaPar_TempOverlap_CU1000_Noint_2cov, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU1000_wide_TempOverlap)
# 
# $sp1_region
# diff         lwr      upr     p adj
# Tropics-N. Temperate      0.8027855 -0.01638814 1.621959 0.0561990
# S. Temperate-N. Temperate 1.0376450  0.50423777 1.571052 0.0000180
# S. Temperate-Tropics      0.2348595 -0.66316475 1.132884 0.8121353


# Non parametric ANCOVA ----




# All angiosperms with Eudicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots$sp1_n_collections + species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots)

Anova(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value  Pr(>F)    
# (Intercept)         12141.7   1 2588.0098 < 2e-16 ***
#   total_n_collections     4.5   1    0.9608 0.32745    
# pair_age                0.6   1    0.1188 0.73050    
# sp1_region             41.8   2    4.4543 0.01208 *  
#   pair_age:sp1_region    47.3   2    5.0458 0.00676 ** 
#   Residuals            2402.1 512                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend     SE  df lower.CL upper.CL
# N. Temperate        -0.1620 0.0389 512  -0.2384  -0.0855
# S. Temperate        -0.0263 0.0326 512  -0.0904   0.0378
# Tropics              0.1410 0.1273 512  -0.1091   0.3911
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate     SE  df t.ratio p.value
# N. Temperate - S. Temperate   -0.136 0.0508 512  -2.669  0.0214
# N. Temperate - Tropics        -0.303 0.1331 512  -2.276  0.0601
# S. Temperate - Tropics        -0.167 0.1314 512  -1.273  0.4109
# 
# P value adjustment: tukey method for comparing a family of 3 estimates  

# NT vs ST is significant and NT vs Tr is almost significant (0.06).

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))

# Because region is significant, let's run a model without interactions considering only regions as main effect
ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots)
Anova(ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections    6.50   1  1.3643   0.24333    
# pair_age              42.86   1  8.9942   0.00284 ** 
#   sp1_region           126.78   2 13.3018 2.332e-06 ***
#   Residuals           2449.41 514                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Tukey with region as main effect
TukeyHSD(ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots)
# 
# $sp1_region
# diff         lwr      upr     p adj
# Tropics-N. Temperate      0.7979930 -0.02260773 1.618594 0.0586972
# S. Temperate-N. Temperate 1.0424754  0.50643730 1.578513 0.0000181
# S. Temperate-Tropics      0.2444824 -0.65590764 1.144872 0.7990833

# Non parametric ANCOVA ----




# Test by taxonomic groups ----

# Monocots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU1000_wide_TempOverlap_monocots$sp1_n_collections + species_pairs_bioclim_CU1000_wide_TempOverlap_monocots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU1000_int_2cov_monocots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_monocots)

Anova(ancovaPar_TempOverlap_CU1000_int_2cov_monocots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq Df  F value  Pr(>F)    
# (Intercept)         2152.46  1 415.8320 < 2e-16 ***
#   total_n_collections    0.00  1   0.0001 0.99320    
# pair_age               0.58  1   0.1118 0.73890    
# sp1_region            11.41  2   1.1025 0.33630    
# pair_age:sp1_region   40.50  2   3.9121 0.02334 *  
#   Residuals            486.57 94                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU1000_int_2cov_monocots, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend     SE df lower.CL upper.CL
# N. Temperate        -0.2727 0.0961 94   -0.463  -0.0819
# S. Temperate         0.3733 0.2162 94   -0.056   0.8027
# Tropics              0.0158 0.2571 94   -0.495   0.5264
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate    SE df t.ratio p.value
# N. Temperate - S. Temperate   -0.646 0.238 94  -2.717  0.0212
# N. Temperate - Tropics        -0.288 0.274 94  -1.051  0.5466
# S. Temperate - Tropics         0.358 0.336 94   1.064  0.5389
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# NT vs ST is significant and NT vs Tr is almost significant (0.06).

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU1000_int_2cov_monocots, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))


# Non parametric ANCOVA ----



# Dicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU1000_wide_TempOverlap_dicots$sp1_n_collections + species_pairs_bioclim_CU1000_wide_TempOverlap_dicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU1000_int_2cov_dicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_dicots)

Anova(ancovaPar_TempOverlap_CU1000_int_2cov_dicots, type = 3 )

# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         9351.3   1 2060.9926 < 2.2e-16 ***
#   total_n_collections    5.4   1    1.1916 0.2756458    
# pair_age               0.0   1    0.0003 0.9867450    
# sp1_region            67.1   2    7.3922 0.0007011 ***
#   pair_age:sp1_region   28.0   2    3.0825 0.0469024 *  
#   Residuals           1873.9 413                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU1000_int_2cov_dicots, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend     SE  df lower.CL upper.CL
# N. Temperate        -0.1395 0.0424 413   -0.223 -0.05623
# S. Temperate        -0.0505 0.0305 413   -0.110  0.00942
# Tropics              0.1874 0.1466 413   -0.101  0.47550
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate     SE  df t.ratio p.value
# N. Temperate - S. Temperate   -0.089 0.0522 413  -1.705  0.2045
# N. Temperate - Tropics        -0.327 0.1526 413  -2.143  0.0826
# S. Temperate - Tropics        -0.238 0.1497 413  -1.589  0.2515
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# None is significant

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU1000_int_2cov_dicots, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))

# Because region was significant, run a model without interactions and then a Tukey
# model without interaction
ancovaPar_TempOverlap_CU1000_Noint_2cov_dicots = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_dicots)
Anova(ancovaPar_TempOverlap_CU1000_Noint_2cov_dicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections    6.37   1  1.3897  0.239125    
# pair_age              41.17   1  8.9826  0.002889 ** 
#   sp1_region           139.01   2 15.1663 4.395e-07 ***
#   Residuals           1901.87 415                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Both age and region are significant. 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU1000_Noint_2cov_dicots, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU1000_wide_TempOverlap_dicots)
# 
# $sp1_region
# diff         lwr      upr     p adj
# Tropics-N. Temperate      0.9727005  0.04628635 1.899115 0.0369996
# S. Temperate-N. Temperate 1.1786307  0.60296852 1.754293 0.0000061
# S. Temperate-Tropics      0.2059302 -0.79896173 1.210822 0.8798571


# Non parametric ANCOVA ----



# Eudicots ----

# Parametric ANCOVA ----
# Get total number of samples
total_n_collections = c(species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots$sp1_n_collections + species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots$sp2_n_collections)

# Set contrasts to be orthogonal: By specifying orthogonal contrasts, we will estimate the the sum of squares for the covariates to be completely partitioned and non-overlapping. In this way we can see the variation attributed to each predictor cleanly and clearly. This is relevant particularly when we have interactions
op = options(contrasts = c("contr.helmert", "contr.poly"))

# Run model
ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots_eudicots = aov(temperature_overlap ~ total_n_collections + pair_age * sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots)

Anova(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots_eudicots, type = 3 )
# Anova Table (Type III tests)
# 
# Response: temperature_overlap
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         9311.4   1 2046.4691 < 2.2e-16 ***
#   total_n_collections    5.3   1    1.1605  0.281986    
# pair_age               0.0   1    0.0014  0.969966    
# sp1_region            63.0   2    6.9246  0.001102 ** 
#   pair_age:sp1_region   29.8   2    3.2781  0.038688 *  
#   Residuals           1870.0 411                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interaction is significant! We need to run a post hoc analysis in which predictors intereact. For this, we will use the library emmeans https://cran.r-project.org/web/packages/emmeans/ Check here https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates

emtrends(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots_eudicots, pairwise ~ sp1_region, var = "pair_age")
# $emtrends
# sp1_region   pair_age.trend     SE  df lower.CL upper.CL
# N. Temperate        -0.1400 0.0424 411   -0.223  -0.0566
# S. Temperate        -0.0415 0.0327 411   -0.106   0.0229
# Tropics              0.1874 0.1468 411   -0.101   0.4759
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate     SE  df t.ratio p.value
# N. Temperate - S. Temperate  -0.0986 0.0536 411  -1.838  0.1586
# N. Temperate - Tropics       -0.3274 0.1528 411  -2.143  0.0826
# S. Temperate - Tropics       -0.2289 0.1504 411  -1.522  0.2817
# 
# P value adjustment: tukey method for comparing a family of 3 estimates  

# None is significant

# Let's visualize slopes
# emms = emmip(ancovaPar_TempOverlap_CU1000_int_2cov_onlyEudicots_eudicots, sp1_region ~ pair_age, cov.reduce = range, plotit = F)
# emmip_ggplot(emms = emms) +theme_classic() + xlab("Divergence") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), name="Region", labels=c("N. Temperate", "S. Temperate", "Tropics"))

# Because region was significant, run a model without interactions and then a Tukey
# model without interaction
ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots_eudicots = aov(temperature_overlap ~ total_n_collections  + pair_age + sp1_region, species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots)
Anova(ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots_eudicots, type = 2)
# Anova Table (Type II tests)
# 
# Response: temperature_overlap
# Sum Sq  Df F value    Pr(>F)    
# total_n_collections    6.38   1  1.3873  0.239546    
# pair_age              34.49   1  7.4973  0.006446 ** 
#   sp1_region           138.29   2 15.0314 4.996e-07 ***
#   Residuals           1899.87 413                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Both age and region are significant. 

# Run post hoc test using region as explanatory variable
TukeyHSD(ancovaPar_TempOverlap_CU1000_Noint_2cov_onlyEudicots_eudicots, "sp1_region", ordered  = TRUE)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# factor levels have been ordered
# 
# Fit: aov(formula = temperature_overlap ~ total_n_collections + pair_age + sp1_region, data = species_pairs_bioclim_CU1000_wide_TempOverlap_onlyEudicots_eudicots)
# 
# $sp1_region
# diff         lwr      upr     p adj
# Tropics-N. Temperate      0.9663469  0.03799071 1.894703 0.0391114
# S. Temperate-N. Temperate 1.1851018  0.60604592 1.764158 0.0000062
# S. Temperate-Tropics      0.2187549 -0.78921275 1.226723 0.8662816



# Non parametric ANCOVA ----








#-----------------------------------------------------------------------#
#-----------------------END OF ANALYSES PREDICTION 2--------------------#
#-----------------------------------------------------------------------#

# Remove objects from workspace to save ====

#---------------------------------#
#                                 #
#            CLEAN UP             #
#                                 #
#---------------------------------#

## Remove unimportant large objects. This helps make the RData image more compact

rm("all_monocots_CU10", 
   "all_monocots_CU100", 
   "all_monocots_CU500", 
   "all_monocots_CU1000", 
   "all_dicots_CU10", 
   "all_dicots_CU100", 
   "all_dicots_CU500", 
   "all_dicots_CU1000",
   "tmp_ntemp_dicots_CU1000_perSp",
   "tmp_ntemp_dicots_CU500_perSp",
   "tmp_ntemp_monocots_CU1000_perSp",
   "tmp_ntemp_monocots_CU500_perSp",
   "tmp_stemp_dicots_CU1000_perSp",
   "tmp_stemp_dicots_CU500_perSp",
   "tmp_stemp_monocots_CU1000_perSp",
   "tmp_stemp_monocots_CU500_perSp",
   "tmp_trop_dicots_CU1000_perSp",
   "tmp_trop_dicots_CU500_perSp",
   "tmp_trop_monocots_CU1000_perSp",
   "tmp_trop_monocots_CU500_perSp",
   "tmp_trop_dicots_CU500_perSp_NoMex",
   "tmp_trop_dicots_CU1000_perSp_NoMex",
   "tmp_trop_monocots_CU500_perSp_NoMex",
   "tmp_trop_monocots_CU1000_perSp_NoMex",
   "tmp_ntemp_monocots_CU500_bioclim_perSp",
   "tmp_trop_monocots_CU500_bioclim_perSp",
   "tmp_stemp_monocots_CU500_bioclim_perSp",
   "tmp_ntemp_dicots_CU500_bioclim_perSp",
   "tmp_trop_dicots_CU500_bioclim_perSp",
   "tmp_stemp_dicots_CU500_bioclim_perSp",
   "tmp_ntemp_monocots_CU1000_bioclim_perSp",
   "tmp_trop_monocots_CU1000_bioclim_perSp",
   "tmp_stemp_monocots_CU1000_bioclim_perSp",
   "tmp_ntemp_dicots_CU1000_bioclim_perSp",
   "tmp_trop_dicots_CU1000_bioclim_perSp",
   "tmp_stemp_dicots_CU1000_bioclim_perSp",
   "genera_ntemp_monocots_CU500",
   "genera_stemp_monocots_CU500",
   "genera_trop_monocots_CU500",
   "genera_ntemp_dicots_CU500",
   "genera_stemp_dicots_CU500",
   "genera_trop_dicots_CU500",
   "genera_ntemp_monocots_CU1000",
   "genera_stemp_monocots_CU1000",
   "genera_trop_monocots_CU1000",
   "genera_ntemp_dicots_CU1000",
   "genera_stemp_dicots_CU1000",
   "genera_trop_dicots_CU1000",
   "all_monocots_CU500_elevRange_coll_group_region",
   "all_dicots_CU500_elevRange_coll_group_region",
   "all_monocots_CU1000_elevRange_coll_group_region",
   "all_dicots_CU1000_elevRange_coll_group_region",
   "all_angiosperms_CU500_elevRange_coll_group_region",
   "all_angiosperms_CU1000_elevRange_coll_group_region",
   "all_monocots_CU500_tempRange_precipRange_coll_group_region",
   "all_dicots_CU500_tempRange_precipRange_coll_group_region",
   "all_monocots_CU1000_tempRange_precipRange_coll_group_region",
   "all_dicots_CU1000_tempRange_precipRange_coll_group_region",
   "all_angiosperms_CU1000_tempRange_precipRange_coll_group_region",
   "all_angiosperms_CU500_tempRange_precipRange_coll_group_region",
   "ntemp_monocots_CU500_bioclim",
   "trop_monocots_CU500_bioclim",
   "stemp_monocots_CU500_bioclim",
   "ntemp_dicots_CU500_bioclim",
   "trop_dicots_CU500_bioclim",
   "stemp_dicots_CU500_bioclim",
   "ntemp_monocots_CU1000_bioclim",
   "trop_monocots_CU1000_bioclim",
   "stemp_monocots_CU1000_bioclim",
   "ntemp_dicots_CU1000_bioclim",
   "trop_dicots_CU1000_bioclim",
   "stemp_dicots_CU1000_bioclim",
   "all_monocots_CU500_NoMex_elevRange_coll_group_region",
   "all_dicots_CU500_NoMex_elevRange_coll_group_region",
   "all_monocots_CU1000_NoMex_elevRange_coll_group_region",
   "all_dicots_CU1000_NoMex_elevRange_coll_group_region",
   "temp_all_species_pairs_CU500_fullGeoData_perSp",
   "temp_all_species_pairs_CU1000_fullGeoData_perSp",
   "nodes_genera_CU500",
   "nodes_genera_CU1000",
   "species_pairs_CU500_not_inData",
   "species_pairs_CU1000_not_inData",
   "columns_gbif",
   "all_genera_CU500",
   "all_genera_CU1000",
   "tmp_ntemp_dicots_CU1000_bioclim_perSp_onlyEudicots",
   "tmp_ntemp_dicots_CU1000_perSp_onlyEudicots",
   "tmp_ntemp_dicots_CU500_bioclim_perSp_onlyEudicots",
   "tmp_ntemp_dicots_CU500_perSp_onlyEudicots",
   "tmp_stemp_dicots_CU1000_bioclim_perSp_onlyEudicots",
   "tmp_stemp_dicots_CU1000_perSp_onlyEudicots",
   "tmp_stemp_dicots_CU500_bioclim_perSp_onlyEudicots",
   "tmp_stemp_dicots_CU500_perSp_onlyEudicots",
   "tmp_trop_dicots_CU1000_bioclim_perSp_onlyEudicots",
   "tmp_trop_dicots_CU1000_perSp_onlyEudicots",
   "tmp_trop_dicots_CU500_bioclim_perSp_onlyEudicots",
   "tmp_trop_dicots_CU500_perSp_onlyEudicots",
   "all_genera_CU1000_onlyEudicots",
   "all_genera_CU500_onlyEudicots",
   "all_species_pairs_CU1000",
   "all_species_pairs_CU1000_fullGeoData",
   "all_species_pairs_CU1000_fullGeoData_onlyEudicots",
   "all_species_pairs_CU1000_onlyEudicots",
   "all_species_pairs_CU1000_overlap",
   "all_species_pairs_CU1000_overlap_long",
   "all_species_pairs_CU1000_overlap_long_onlyEudicots",
   "all_species_pairs_CU1000_overlap_onlyEudicots",
   "all_species_pairs_CU500",
   "all_species_pairs_CU500_fullGeoData",
   "all_species_pairs_CU500_fullGeoData_onlyEudicots",
   "all_species_pairs_CU500_onlyEudicots",
   "all_species_pairs_CU500_overlap",
   "all_species_pairs_CU500_overlap_long",
   "all_species_pairs_CU500_overlap_long_onlyEudicots",
   "all_species_pairs_CU500_overlap_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU1000",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData_counts_means_global",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData_counts_means_global_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU500_long_allData_perGroup_perRegion",
   "all_species_pairs_match_all_angiosperms_CU500_long_allData_perGroup_perRegion_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU500_onlyEudicots",
   "all_angiosperms_CU1000_allData_overlap",
   "all_angiosperms_CU1000_allData_overlap_onlyEudicots",
   "all_angiosperms_CU1000_elevRange_coll_group_region_onlyEudicots",
   "all_angiosperms_CU1000_tempRange_precipRange_coll_group_region_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData_perGroup_perRegion",
   "all_species_pairs_match_all_angiosperms_CU1000_long_allData_perGroup_perRegion_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU1000_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU500",                                        "all_species_pairs_match_all_angiosperms_CU500_long_allData",                          "all_species_pairs_match_all_angiosperms_CU500_long_allData_counts_means_global",
   "all_species_pairs_match_all_angiosperms_CU500_long_allData_counts_means_global_onlyEudicots",
   "all_species_pairs_match_all_angiosperms_CU500_long_allData_onlyEudicots",
   "genera_ntemp_dicots_CU1000_onlyEudicots",
   "genera_ntemp_dicots_CU500_onlyEudicots",
   "genera_stemp_dicots_CU1000_onlyEudicots",
   "genera_stemp_dicots_CU500_onlyEudicots",
   "genera_trop_dicots_CU1000_onlyEudicots",
   "genera_trop_dicots_CU500_onlyEudicots",
   "nodes_genera_CU1000_onlyEudicots",
   "nodes_genera_CU500_onlyEudicots",
   "nodes.info.2_tbl",
   "non_eudicots",
   "ntemp_dicots_CU1000_bioclim_onlyEudicots",
   "ntemp_dicots_CU500_bioclim_onlyEudicots",
   "species_pairs_bioclim_CU1000",
   "species_pairs_bioclim_CU1000_onlyEudicots",
   "species_pairs_bioclim_CU1000_wide",
   "species_pairs_bioclim_CU1000_wide_onlyEudicots",
   "species_pairs_bioclim_CU500",
   "species_pairs_bioclim_CU500_onlyEudicots",
   "species_pairs_bioclim_CU500_wide",
   "species_pairs_bioclim_CU500_wide_onlyEudicots",
   "species_pairs_CU1000_not_inData_onlyEudicots",
   "species_pairs_CU500_not_inData_onlyEudicots",
   "species_pairs_elevation_range_CU1000",
   "species_pairs_elevation_range_CU1000_onlyEudicots",
   "species_pairs_elevation_range_CU1000_wide",
   "species_pairs_elevation_range_CU1000_wide_onlyEudicots",
   "species_pairs_elevation_range_CU500",
   "species_pairs_elevation_range_CU500_onlyEudicots",
   "species_pairs_elevation_range_CU500_wide",
   "species_pairs_elevation_range_CU500_wide_onlyEudicots",
   "stemp_dicots_CU1000_bioclim_onlyEudicots",
   "stemp_dicots_CU500_bioclim_onlyEudicots",
   "temp_all_species_pairs_CU1000_fullGeoData_perSp_onlyEudicots",
   "temp_all_species_pairs_CU500_fullGeoData_perSp_onlyEudicots",
   "trop_dicots_CU1000_bioclim_onlyEudicots",
   "trop_dicots_CU500_bioclim_onlyEudicots",
   "species_pairs_elevation_range_CU500_wide_ElevOverlap_monocots",
   "species_pairs_elevation_range_CU500_wide_ElevOverlap_dicots",
   "species_pairs_bioclim_CU500_wide_TempOverlap_monocots",
   "species_pairs_bioclim_CU500_wide_TempOverlap_dicots",
   "species_pairs_elevation_range_CU1000_wide_ElevOverlap_monocots",
   "species_pairs_elevation_range_CU1000_wide_ElevOverlap_dicots",
   "species_pairs_bioclim_CU1000_wide_TempOverlap_monocots",
   "species_pairs_bioclim_CU1000_wide_TempOverlap_dicots",
)


# Record information about the session
session_info_kernel = sessionInfo()
system_time_kernel = Sys.time()

time_stop = Sys.time()
time_run = time_stop - time_start

time_run

# Write the results to prepare them for manuscript.rmd

save.image("manuscript.RData")
