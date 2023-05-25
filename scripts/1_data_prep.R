## Trait-range manuscript
## Mariana Garcia Criado
## March 2020
## Script 1. Data prep (traits & ranges)


# The input files for this script are not available in the repo given that they include raw data.
# This is a data prep script and is made available for reproducibility and transparency purposes. 
# The final file produced by this script (2022_master_new_superfinal.RData) will be used as input data 
# in the rest of scripts.


## LIBRARIES ----
library(dplyr)
library(cowplot)


## DATA PREP ----

# TRY 3.0 database
# tryv3 <- read.csv("data/trait_clean.csv") # 16285 observations
# [this file is not available in the repo as it contains raw data]
# unique(tryv3$ValueKindName) # Single, Mean, Median, Individual Mean, Plot mean, Site specific mean, na (7 types)


## New cleaned TRY 5.0 database
load("data/mTTTT_clean.RData") # 20032 obs [this file is not available in the repo as it contains raw data]
# unique(mTTTT.clean$AccSpeciesName) # 105 species
# unique(mTTTT.clean$ValueKindName) # Single, Mean, Best estimate, Median, Individual Mean, Plot mean, NA, Site Specific mean



# checking the NAs
summaryNA <- mTTTT.clean %>% 
  group_by(ValueKindName, Dataset) %>% 
  mutate(nobs = length(StdValue))

summNA <- count(summaryNA, ValueKindName, Dataset) # 993 observations in NA, we want to retain those

# Retain only Single and NA (as we know the NAs are single values)
single <- mTTTT.clean %>% filter(ValueKindName %in% c("Single", "Individual Mean") | is.na(ValueKindName)) #18804 



## DATASETS FROM TRY 3.0 ----
#load("data/master_full_new.RData") #object master.full.new [this file is not available in the repo as it contains raw data]
#datasets.new <- as.data.frame(unique(master.full.new$Dataset))

# we are missing ArtDeco, BIOPOP and Causasus databases in TRY 5.0
try.ttt <- load("data/try_ttt_clean.RData") # [this file is not available in the repo as it contains raw data]
datasets.raw.v3 <- as.data.frame(unique(try.ttt.clean$Dataset)) # all 3 missing datasets are in TRY 3.0

# filter these 3 datasets
extra.db <- filter(try.ttt.clean, Dataset == "Causasus Plant Traits Database" | 
                     Dataset == "ArtDeco Database" | 
                     Dataset == "BIOPOP: Functional Traits for Nature Conservation")

# run the vector with our species
name.vector.try5.full <- read.csv("data/name_vector_try5_full.csv") # [this file is not available in the repo as it contains raw data]
name.vector.vector.try5.full <- unique(name.vector.try5.full$unique.master.full.new.sp.) #62 species

# extract data for our species and traits and standardise columns
extra.sp <- extra.db %>% 
  filter(AccSpeciesName %in% name.vector.vector.try5.full) %>% 
  filter(TraitShort == "PlantHeight" | TraitShort == "SeedMass" | TraitShort == "SLA") %>%
  filter(ValueKindName %in% c("Single", "Individual Means") | is.na(ValueKindName)) %>% 
  mutate(ObservationID = NA, ObsDataID = NA, Treatment = NA) %>%
  dplyr::select(Dataset, AccSpeciesName, ObservationID, ObsDataID, ValueKindName, StdValue, 
         Lat, Lon, DataContributor, TraitShort, Source, Treatment, Genus)

# bind the two datasets + additional data (extracted from other non-TRY/TTT sources)
vac.oxyc <- read.csv("data/vac_oxyc_add_height.csv") # [this file is not available in the repo as it contains raw data]
sal.arc <- read.csv("data/salix_arctica_add_seed.csv") # [this file is not available in the repo as it contains raw data]

mTTTT.clean2 <- bind_rows(single, extra.sp, vac.oxyc, sal.arc) #18832 obs


# Check if there is any way of extracting the location of species that don't have coordinate data
na.lat <- mTTTT.clean2 %>% tidyr::drop_na(StdValue) %>% filter(is.na(Lat))
na.lat.group <- na.lat %>% group_by(Dataset) %>% mutate(n())
summ.na.lat <- count(na.lat.group, Dataset)

# I've checked these databases and retained those whose names indicated the location of the records
# Only retained those that clearly were above 30 degrees North

# create vector with database names
na.dbs <- c("Abisko & Sheffield Database", "Causasus Plant Traits Database", "Ecological Flora of the British Isles",
            "Italian Alps Plant Traits Database", "Plant Traits from Romania", "Sheffield Database", 
            "The LEDA Traitbase", "The Netherlands Plant Height Database", "TOPIC (Traits of Plants in Canada)",
            "Traits for Herbaceous Species from Andorra", "Tundra Plant Traits Database", "Urals")

# Filter the coord NA from these datasets
na.chosen <- mTTTT.clean2 %>% filter(is.na(Lat)) %>% filter(Dataset %in% na.dbs)

# prepare extra seed datasets here to bind below
seed.leda <- read.csv("data/seed_mass_extra_LEDA.csv")
seed.lucie <- read.csv("data/seed_mass_lucie_smrzova.csv")
seed.bob <- read.csv("data/seed_mass_hollister.csv")



## FILTER TRAIT DATA ----

# filter for species north of 30 degrees latitude + add extra seed data
try50 <- mTTTT.clean2 %>% 
  tidyr::drop_na(Lat) %>% 
  filter(Lat > 30) %>% 
  bind_rows(na.chosen, seed.leda, seed.lucie, seed.bob) %>%
  rename(sp = AccSpeciesName) #18220

# Single, Individual Mean and NA (which are also single)
unique(try50$ValueKindName) 


## Remove here the problematic measurements so the observation number is accurate below
## It was decided with a few co-authors that these measures are likely to be erroneous and measuring length rather than height
try5 <- try50 %>% 
  filter(!(TraitShort == "PlantHeight" & sp == "Linnaea borealis" & StdValue %in% c(1.2, 0.525, 0.4))) %>%
  filter(!(TraitShort == "PlantHeight" & sp == "Vaccinium oxycoccos" & StdValue %in% c(0.8, 0.6, 0.45, 0.3)))
# 18213 obs



## REMOVE DUPLICATES ----
twice5 <- try5 %>% 
  group_by(sp, TraitShort, Lat, Lon, StdValue) %>% 
  mutate(dupe = n()>1)

# filter per Source
trues.TTT <- twice5 %>% filter(Source == "TTT") %>% filter(dupe == "TRUE") #6821 obs 
trues.TRY <- twice5 %>% filter(Source == "TRY") %>% filter(dupe == "TRUE") #1253 obs

# Are these actual duplicates in TTT
summaryTTT <- trues.TTT %>% 
  group_by(TraitShort, sp, Dataset, DataContributor) %>% 
  mutate(nobs = length(StdValue))

summTTT <- count(summaryTTT, Dataset)


# Are these actual duplicates in TRY
summaryTRY <- trues.TRY %>% 
  group_by(TraitShort, sp, Dataset, DataContributor) %>% 
  mutate(nobs = length(StdValue))

summTRY <- count(summaryTRY, Dataset)

# there are some values of nobs = 1, this is probably because the duplicates are not within TRY but with TTT, 
# which is why I have also tried to identify the duplicates within databases without TRY/TTT distinction


# check duplicates without Source distinction
summaryALL <- twice5 %>% filter(dupe == "TRUE") %>% 
  group_by(TraitShort, Dataset, DataContributor) %>% 
  mutate(nobs = length(StdValue))

summALL <- count(summaryALL, Dataset) 


# Check if there are duplicates included both in TRY and TTT
contribs.ttt <- unique(summTTT$DataContributor) #35 contributors
contribs.try <- unique(summTRY$DataContributor) #29 contributors 

same.contrib <- contribs.try %in% contribs.ttt
# Logan Berner is in both datasets 

# every single value without exception appears both in TTT and TRY
logan <- filter(try5, DataContributor == "Logan Berner")

# 42 values including duplicates, which means half (the ones in TRY) should go (21)
# this is for Betula nana which has almost 700 observations so removing these doesn't affect the min #obs

# Remove the TRY duplicates
try5.clean <- try5 %>% 
  filter(!(DataContributor == "Logan Berner" & Source == "TRY"))
# 18192 obs

# Alastair Rogers' data has duplicates too - leave only unique values (unique ObservationID)
rogers <- filter(try5, DataContributor == "Alistair Rogers") %>% 
  distinct(ObservationID, .keep_all = TRUE)
# that leaves 108 obs - before distinct there were 216

# Salix rotundifolia has 20 observations (that leaves 10 so still works), Salix pulchra has over 200
try5.clean2 <- try5.clean %>% filter(!(DataContributor == "Alistair Rogers"))

# bind these two with the clean data
try5.clean3 <- bind_rows(try5.clean2, rogers) #18084 obs



## OBSERVATIONS NUMBER FILTER ----

# Filter those with more than 4 records (min 5) per species and trait  
try5.clean3obs <- try5.clean3 %>% 
  group_by(sp, TraitShort) %>% 
  mutate(nobs = length(StdValue)) %>% 
  filter(nobs > 4) %>%
  ungroup()
## 17940 observations



## FIX ERRORS ----

# there are 6 records from Plant Traits in Pollution Gradients Database where the Plant Height values 
# seem wrong - from 4.1000 to 55.0000. The unit is supposed to be meters, and the highest value is 2.8m, 
# so I think these are supposed to be cm instead of m. I have check this with Madhur Anand and she confirms
# that the values are provided in cm and need to be converted to m.

# extract and transform Madhur's data to cm
vac.cm2 <- try5.clean3obs %>% 
  filter(DataContributor == "Madhur Anand") %>% 
  filter(TraitShort == "PlantHeight") %>%
  mutate(StdValue = StdValue/100)

# remove from the original database the erroneous data
master.full.novac2 <- try5.clean3obs %>% 
  filter(!(DataContributor == "Madhur Anand" & TraitShort == "PlantHeight"))

# bind these two in a proper updated database
try5.clean4 <- bind_rows(master.full.novac2, vac.cm2) #17940


# Another issue - a record from Fritz Schweingruber (55.00000, -24.000000) on Harrimanella hypnoides
# in which plantheight = 0.5m. Those coordinates are in the middle of the Atlantic Ocean and there is no land there.
# I have tracked down Fritz's original database where this record has indeed those coordinates but it is
# associated to Scorsbysund, Greenland. The coordinates for this place are 70, -24 instead.
# Sadly Fritz passed away in January 2020 so I can't double-check these records with him, 
# but these are certainly the right coordinates so changing the Latitude here below. 

# extract wrong value
fritz <- try5.clean4 %>%
  filter(DataContributor == "Fritz Schweingruber") %>% 
  filter(TraitShort == "PlantHeight") %>%
  filter(sp == "Harrimanella hypnoides") %>%
  mutate(Lat = 55.00000) %>%
  mutate(Lat = 70.00000)

# remove from the original database the erroneous data
master.new.nofritz <- try5.clean4 %>% 
  filter(!(DataContributor == "Fritz Schweingruber" & TraitShort == "PlantHeight" & sp == "Harrimanella hypnoides"))

# bind these two in a proper updated database
try5.clean5 <- bind_rows(master.new.nofritz, fritz) #17940 obs


# change Ledum palustre to Rhododendron tomentosum here
try5.clean5$sp[try5.clean5$sp == "Ledum palustre"] <- "Rhododendron tomentosum"
try5.clean5$Genus[try5.clean5$sp == "Rhododendron tomentosum"] <- "Rhododendron"

# save this file - 17940 obs
save(try5.clean5, file = "data/2022_try5_clean5.RData") # [this file is not available in the repo as it contains raw data]




## RANGE DATA ----
range.data <- read.csv("data/SRC_baseline_tabs_2017-04-26_fullrange.csv") # [this file is not available in the repo as it contains raw data]

unique(range.data$sp) #132 species

# Convert species names with dot to spaces to match TRY
range.data$sp <- gsub(".", " ", range.data$sp, fixed = TRUE)

# convert the two hiphenated species
range.data$sp[range.data$sp == "Arctostaphylos uvaursi"] <- "Arctostaphylos uva-ursi"
range.data$sp[range.data$sp == "Vaccinium vitisidaea"] <- "Vaccinium vitis-idaea"

# Species name changes: convert species that have PAF nomenclature to The Plant List nomenclature
# I have done this manually comparing Word docs 
range.data$sp[range.data$sp == "Betula fruticosa"] <- "Betula humilis"
range.data$sp[range.data$sp == "Chamaepericlymenum canadense"] <- "Cornus canadensis"
range.data$sp[range.data$sp == "Swida sericea"] <- "Cornus sericea"
range.data$sp[range.data$sp == "Kalmia procumbens"] <- "Loiseleuria procumbens"
range.data$sp[range.data$sp == "Oxycoccus palustris"] <- "Vaccinium oxycoccos"


# Filter right parameters and convert pixel number into km (each pixel size is 10x10km)
# Divide range changes by a number so the values are smaller and log-transforming works properly
range.full <- range.data %>% 
  filter(area == "full_area") %>% 
  filter(biointer == "no") %>%
  filter(filt == "max_dipersal") %>%
  filter(sp %in% name.vector.try5) %>%
  mutate(scenario = paste0(rcp, "_", gcm)) %>%
  mutate(RangeLog = log(CurrentRangeSize)) %>%
  mutate(CurrentRangeSizeKm = CurrentRangeSize*100) %>%
  mutate(RangeLogKm = log(CurrentRangeSizeKm)) %>%
  mutate(FutureRangeSize.FullDisp.Km = FutureRangeSize.FullDisp*100) %>% 
  mutate(AbsoluteRangeChangeKm = FutureRangeSize.FullDisp.Km - CurrentRangeSizeKm) %>%
  mutate(AbsoluteRangeChangeKmDivided = AbsoluteRangeChangeKm/1000000) %>%
  mutate(SpeciesRangeChangeDivided = SpeciesRangeChange/100)

# save this for script#3
write.csv(range.full, file = "data/2022_range_full.csv") # [this file is not available in the repo as it contains raw data]


# Find the smallest value so we can use it to add as a constant
min.abs <- range.full %>% dplyr::select(AbsoluteRangeChangeKmDivided) %>% min() %>% abs()
min.rel <- range.full %>% dplyr::select(SpeciesRangeChangeDivided) %>% min() %>% abs()
min.test <- range.full %>% dplyr::select(SpeciesRangeChange) %>% min() %>% abs()
min.abs.test <- range.full %>% dplyr::select(AbsoluteRangeChangeKm) %>% min() %>% abs()

# Add constant so we can log-transform negative values, calculate quantiles, log-transform and centre data
range.full2 <- range.full %>%
  mutate(AbsoluteRangeChangeKmConstant = AbsoluteRangeChangeKmDivided + min.abs + 1) %>%
  mutate(SpeciesRangeChangeConstant = SpeciesRangeChangeDivided + min.rel + 1) %>%
  group_by(sp) %>% 
  mutate(rel.quan25 = quantile(SpeciesRangeChangeConstant, probs = 0.25)) %>%
  mutate(rel.median = quantile(SpeciesRangeChangeConstant, probs = 0.50)) %>%
  mutate(rel.quan75 = quantile(SpeciesRangeChangeConstant, probs = 0.75)) %>%
  mutate(abs.quan25 = quantile(AbsoluteRangeChangeKmConstant, probs = 0.25)) %>%
  mutate(abs.median = quantile(AbsoluteRangeChangeKmConstant, probs = 0.50)) %>%
  mutate(abs.quan75 = quantile(AbsoluteRangeChangeKmConstant, probs = 0.75)) %>%
  mutate(future.med = quantile(FutureRangeSize.FullDisp.Km, probs = 0.50)) %>%
  ungroup() %>%
  mutate(rel.quan25.log = log(rel.quan25)) %>%
  mutate(rel.median.log = log(rel.median)) %>%
  mutate(rel.quan75.log = log(rel.quan75)) %>%
  mutate(abs.quan25.log = log(abs.quan25)) %>%
  mutate(abs.median.log = log(abs.median)) %>%
  mutate(abs.quan75.log = log(abs.quan75)) %>%
  distinct(sp, .keep_all = TRUE) %>%
  mutate(CentredRangeLog = RangeLogKm - mean(RangeLogKm)) %>%
  mutate(cent.rel.quan25.log = rel.quan25.log - mean(rel.quan25.log)) %>%
  mutate(cent.rel.median.log = rel.median.log - mean(rel.median.log)) %>%
  mutate(cent.rel.quan75.log = rel.quan75.log - mean(rel.quan75.log)) %>%
  mutate(cent.abs.quan25.log = abs.quan25.log - mean(abs.quan25.log)) %>%
  mutate(cent.abs.median.log = abs.median.log - mean(abs.median.log)) %>%
  mutate(cent.abs.quan75.log = abs.quan75.log - mean(abs.quan75.log))

# Species with trait data - 62 observations, one per species with the clim scenarios condensed on quantiles.

# At some point I got a warning saying NaNs have been created so just confirming this is not the case:
# any(is.nan(range.full$rel.quan25.log))
# any(is.nan(range.full$rel.median.log))
# any(is.nan(range.full$rel.quan75.log))
# any(is.nan(range.full$abs.quan25.log))
# any(is.nan(range.full$abs.median.log))
# any(is.nan(range.full$abs.quan75.log))


# Checking distributions
hist(range.full2$cent.rel.quan25.log)
hist(range.full2$rel.quan25.log)
hist(range.full2$cent.rel.median.log)
hist(range.full2$rel.median.log)
hist(range.full2$cent.rel.quan75.log)
hist(range.full2$rel.quan75.log)

hist(range.full2$cent.abs.quan25.log)
hist(range.full2$abs.quan25.log)
hist(range.full2$cent.abs.median.log)
hist(range.full2$abs.median.log)
hist(range.full2$cent.abs.quan75.log)
hist(range.full2$abs.quan75.log)

# save this for later
write.csv(range.full2, file = "data/2022_range_full2.csv")
# this file is saved with all the columns but I have kept only those relevant to the analyses and removed the raw data


## MERGING DATABASES ----
master.data <- merge(try5.clean5, range.full2, by = "sp")
# 17921 obs

vector.merge <- unique(master.data$sp) # 62 species
name.vector.try5.2 <- unique(try5.clean5$sp) # 63 species

similar.spps <- vector.merge %in% name.vector.try5
similar.spps2 <- name.vector.try5 %in% vector.merge

## Ribes rubrum is the species in the original try5 vector that does not have range data
## Ribes rubrum does not appear as a species on its own on the range data, instead we have Ribes spicatum,
## and the PAF mastersheets consider R. rubrum to be a synonym of R. spicatum. 
## However TRY considers these to be 2 different species (and so does the Plant List)
## and there are TRY data for both. Since it is not possible to fully clarify if the range data refers to R. rubrum
## we are going to have to lose the trait data for R. rubrum and just go with R. spicatum.

# check if this accounts for the amount of records lost
ribs <- filter(try5.clean5, sp == "Ribes rubrum")
# checks out - these are the 19 records lost in the merge - all good!



## ADD RANGE CATEGORY ----

# Adding qualitative range categories
range.cat.hist <- master.data %>% distinct(sp, .keep_all = TRUE)

# histogram of current ranges
hist(range.cat.hist$CurrentRangeSizeKm, breaks = 100)

# calculate quantiles to define categories
quantile(range.cat.hist$CurrentRangeSizeKm, 0.2) #3092760
quantile(range.cat.hist$CurrentRangeSizeKm, 0.5) #6909000
quantile(range.cat.hist$CurrentRangeSizeKm, 0.8) #10686760

# create another object
range.trait <- master.data

# add qualitative categories as a function of current range
range.trait$range_cat[range.trait$CurrentRangeSizeKm < 3092760] <- "small"
range.trait$range_cat[range.trait$CurrentRangeSizeKm >= 3092760 & 
                        range.trait$CurrentRangeSizeKm <= 10686760] <- "medium"
range.trait$range_cat[range.trait$CurrentRangeSizeKm > 10686760] <- "large"

# without duplicates: L-13, M-36, S-13
cat.count <- distinct(range.trait, sp, .keep_all = TRUE) %>% count(range_cat)




## LEAF PHENOLOGY ----
form <- read.csv("data/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.csv") 
# [this file is not available in the repo as it contains raw data]

# filter for the relevant species in our list
form.data <- form %>% filter(AccSpeciesName %in% vector.merge) 
form.short <- dplyr::select(form.data, AccSpeciesName, Family, LeafPhenology)
form.short$AccSpeciesName # data for 54 species

# join with main trait database
names(form.short)[names(form.short) == "AccSpeciesName"] <- "sp"
master.full <- left_join(range.trait, form.short, by = "sp")

# which species have missing data?
unique(master.full$LeafPhenology)
leaf.na <- subset(master.full, is.na(LeafPhenology))
unique(leaf.na$sp)

# empty rows
leaf.empty <- subset(master.full, LeafPhenology == "")
unique(leaf.empty$sp)

# gap-fill the missing data for these species
master.full$LeafPhenology[master.full$sp == "Arctous rubra"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Dasiphora fruticosa"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Harrimanella hypnoides"] <- "evergreen"
master.full$LeafPhenology[master.full$sp == "Salix argyrocarpa"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Salix vestita"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Rhododendron tomentosum"] <- "evergreen"

master.full$LeafPhenology[master.full$sp == "Artemisia dracunculus"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Artemisia gmelinii"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Ribes glandulosum"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Salix barrattiana"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Salix niphoclada"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Salix phlebophylla"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Salix rotundifolia"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Sibbaldia procumbens"] <- "evergreen"
master.full$LeafPhenology[master.full$sp == "Vaccinium caespitosum"] <- "deciduous"
master.full$LeafPhenology[master.full$sp == "Vaccinium myrtilloides"] <- "deciduous"

# M. gale was coded as evergreen/deciduous - seems to be deciduous
master.full$LeafPhenology[master.full$sp == "Myrica gale"] <- "deciduous"

# How many levels per leaf phenology?
leaf.group <- master.full %>% group_by(LeafPhenology) %>% count()
(unique(master.full$LeafPhenology))

# Modify with capitals for the graphs
master.full$LeafPhenology <- as.character(master.full$LeafPhenology)

master.full$LeafPhenology[master.full$LeafPhenology == "evergreen"] <- "Evergreen"
master.full$LeafPhenology[master.full$LeafPhenology == "deciduous"] <- "Deciduous"
(unique(master.full$LeafPhenology))



## FAMILY DATA ----

# which sp have missing family data?
unique(master.full$Family)
fam.na <- subset(master.full, is.na(Family))
unique(fam.na$sp) # 8 species

# empty rows
fam.empty <- subset(master.full, Family == "")
unique(fam.empty$sp) # none

# gap-fill families
master.full$Family[master.full$sp == "Arctous rubra"] <- "Ericaceae"
master.full$Family[master.full$sp == "Dasiphora fruticosa"] <- "Rosaceae"
master.full$Family[master.full$sp == "Harrimanella hypnoides"] <- "Ericaceae"
master.full$Family[master.full$sp == "Vaccinium myrtilloides"] <- "Ericaceae"
master.full$Family[master.full$sp == "Salix vestita"] <- "Salicaceae"
master.full$Family[master.full$sp == "Salix argyrocarpa"] <- "Salicaceae"
master.full$Family[master.full$sp == "Ribes glandulosum"] <- "Grossulariaceae"
master.full$Family[master.full$sp == "Rhododendron tomentosum"] <- "Ericaceae"
           
# How many levels per family?
fam.group <- master.full %>% group_by(Family) %>% mutate(unique.sp = n_distinct(sp))
fam.sum <- fam.group %>% distinct(Family, .keep_all = TRUE) %>% dplyr::select(Family, unique.sp)
# 3 families out of 12 have more than 3 levels (species). 3 families have 3 levels.
# 6 families have less than 3 levels (species).



## FUNCTIONAL GROUPS ----
# Our FG data was coded by Anne BO and Signe as follows: "For each species, we computed the median of the maximum 
# height (in meters) to classify the shrubs into different height classes following the classification 
# of Myers-Smith et al. (2015a): dwarf shrubs (< 0.2 m), low shrubs (0.2-0.5 m), and tall shrubs (> 0.5 m)" with data
# from TRY, but some of these categories are off. So we have cross-checked all the maximum heights with online flora
# (see cat_traits_isla_IMS_MGC.xlsx and Shrub_canopy_heights_funcgroup_MGC) and I am re-coding the categories here.

# Modifying some functional groups that seem to be erroneous
master.full$gf[master.full$sp == "Andromeda glaucophylla"] <- "Low shrub"
master.full$gf[master.full$sp == "Andromeda polifolia"] <- "Tall shrub"
master.full$gf[master.full$sp == "Arctostaphylos uva-ursi"] <- "Dwarf shrub"
master.full$gf[master.full$sp == "Arctous alpina"] <- "Low shrub"
master.full$gf[master.full$sp == "Arctous rubra"] <- "Dwarf shrub"
master.full$gf[master.full$sp == "Artemisia frigida"] <- "Low shrub"
master.full$gf[master.full$sp == "Betula nana"] <- "Tall shrub"
master.full$gf[master.full$sp == "Cassiope tetragona"] <- "Low shrub"
master.full$gf[master.full$sp == "Empetrum nigrum"] <- "Low shrub"
master.full$gf[master.full$sp == "Loiseleuria procumbens"] <- "Low shrub"
master.full$gf[master.full$sp == "Ribes glandulosum"] <- "Low shrub"
master.full$gf[master.full$sp == "Rubus chamaemorus"] <- "Low shrub"
master.full$gf[master.full$sp == "Salix arctica"] <- "Low shrub"
master.full$gf[master.full$sp == "Salix lanata"] <- "Tall shrub"
master.full$gf[master.full$sp == "Salix myrsinites"] <- "Tall shrub"
master.full$gf[master.full$sp == "Vaccinium caespitosum"] <- "Tall shrub"
master.full$gf[master.full$sp == "Vaccinium myrtillus"] <- "Tall shrub"
master.full$gf[master.full$sp == "Vaccinium uliginosum"] <- "Tall shrub"




## DISPERSAL MODE ----

# Check our species again
unique(master.full$sp)

# Collect berry species in a vector
berry.vector <- c("Vaccinium caespitosum", "Vaccinium myrtilloides", "Vaccinium myrtillus", 
                  "Vaccinium uliginosum", "Vaccinium vitis-idaea", "Vaccinium oxycoccos", 
                  "Empetrum nigrum", "Rubus chamaemorus",
                  "Rubus idaeus", "Cornus canadensis", "Arctostaphylos uva-ursi", "Arctous alpina",
                  "Arctous rubra", "Ribes spicatum", "Ribes glandulosum", "Ribes nigrum")

# Add in dispersal mode per species
master.new.superfinal <- mutate(master.full, DispersalMode = ifelse(sp %in% berry.vector, "Berry", "Wind"))

# Double-check that they are all coded
unique(master.new.superfinal$DispersalMode) # all ok

# keep only relevant columns
master.new.superfinal <- master.new.superfinal %>% 
  select(-c(X, layer_id, Loss, Stable0, Stable1, Gain, PercLoss, PercGain, CurrentRangeSize,
            FutureRangeSize.NoDisp, FutureRangeSize.FullDisp, area, file.id, model, rcp, gcm, biointer, src_ras_file, 
            scenario, RangeLog, CurrentRangeSizeKm, RangeLogKm, FutureRangeSize.FullDisp.Km, AbsoluteRangeChangeKmDivided, 
            SpeciesRangeChangeDivided, AbsoluteRangeChangeKmConstant, SpeciesRangeChangeConstant))

# save database
save(master.new.superfinal, file = "data/2022_master_new_superfinal.RData")
# 17921 obs


# BEFORE MOVING ON...
#each row should always have a different combination of trait, trait value or latitude/longitude
#for each species, the range values should be the same
test.harry <- filter(master.new.superfinal, sp == "Harrimanella hypnoides") # 107 obs





## ITV METHOD COMPARISON ----

# Compare the different trait values and variation from usign different trait metrics
# These become Supplementary Data files 2, 3 and 4.

load("data/2022_master_new_superfinal.RData")

# prepare dataset
itv.db <- master.new.superfinal %>% select(sp, Dataset, StdValue, Lat, Lon, TraitShort) %>%
  tidyr::unite('Location', Lat:Lon, remove = F, sep = ";") %>% 
  mutate(StdValue = case_when(TraitShort == "SeedMass" & StdValue == 0.000000000 ~ 0.000000001, TRUE ~ StdValue)) %>%
  mutate(LogTraitValue = log(StdValue))


# SLA with all records
sla.itv <- itv.db %>% filter(TraitShort == "SLA") %>% 
  group_by(sp) %>% 
  mutate(NumbRecords = length(StdValue), NumbSites = length(unique(Location))) %>%
  mutate(MTV = median(StdValue), LogMTV = median(LogTraitValue), 
         Mean = mean(StdValue), LogMean = mean(LogTraitValue),
         SD = sd(StdValue), LogSD = sd(LogTraitValue), 
         COV = SD/Mean, LogCOV = LogSD/LogMean, 
         LogCOVnew = log(SD/Mean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, TraitShort, NumbRecords, NumbSites, MTV, LogMTV, SD, LogSD, COV, LogCOV, LogCOVnew)

# SLA with 5 random records
sla.itv.random <- itv.db %>% filter(TraitShort == "SLA") %>% 
  group_by(sp) %>% 
  slice_sample(n = 5) %>%
  mutate(SampleNumbRecords = 5, SampleNumbSites = length(unique(Location))) %>%
  mutate(SampleMTV = median(StdValue), LogSampleMTV = median(LogTraitValue), 
         SampleMean = mean(StdValue), SampleLogMean = mean(LogTraitValue),
         SampleSD = sd(StdValue), LogSampleSD = sd(LogTraitValue), 
         SampleCOV = SampleSD/SampleMean, SampleLogCOV = LogSampleSD/SampleLogMean, 
         SampleLogCOVnew = log(SampleSD/SampleMean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, SampleNumbRecords, SampleNumbSites, SampleMTV, LogSampleMTV, 
         SampleSD, LogSampleSD, SampleCOV, SampleLogCOV, SampleLogCOVnew)

# Full table
sla.itv.all <- left_join(sla.itv, sla.itv.random, by = "sp")

write.csv(sla.itv.all, "data/2022_sla_itv_comparison.csv")




# Seed mass with all records
seed.itv <- itv.db %>% filter(TraitShort == "SeedMass") %>% 
  group_by(sp) %>% 
  mutate(NumbRecords = length(StdValue), NumbSites = length(unique(Location))) %>%
  mutate(MTV = median(StdValue), LogMTV = median(LogTraitValue), 
         Mean = mean(StdValue), LogMean = mean(LogTraitValue),
         SD = sd(StdValue), LogSD = sd(LogTraitValue), 
         COV = SD/Mean, LogCOV = LogSD/LogMean,
         LogCOVnew = log(SD/Mean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, TraitShort, NumbRecords, NumbSites, MTV, LogMTV, SD, LogSD, COV, LogCOV, LogCOVnew)

# Seed Mass with 5 random records
seed.itv.random <- itv.db %>% filter(TraitShort == "SeedMass") %>% 
  group_by(sp) %>% 
  slice_sample(n = 5) %>%
  mutate(SampleNumbRecords = 5, SampleNumbSites = length(unique(Location))) %>%
  mutate(SampleMTV = median(StdValue), LogSampleMTV = median(LogTraitValue), 
         SampleMean = mean(StdValue), SampleLogMean = mean(LogTraitValue),
         SampleSD = sd(StdValue), LogSampleSD = sd(LogTraitValue), 
         SampleCOV = SampleSD/SampleMean, SampleLogCOV = LogSampleSD/SampleLogMean,
         SampleLogCOVnew = log(SampleSD/SampleMean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, SampleNumbRecords, SampleNumbSites, SampleMTV, LogSampleMTV, 
         SampleSD, LogSampleSD, SampleCOV, SampleLogCOV, SampleLogCOVnew)

# Full table
seed.itv.all <- left_join(seed.itv, seed.itv.random, by = "sp")

write.csv(seed.itv.all, "data/2022_seed_itv_comparison.csv")



# Height with all records
hei.itv <- itv.db %>% filter(TraitShort == "PlantHeight") %>% 
  group_by(sp) %>% 
  mutate(NumbRecords = length(StdValue), NumbSites = length(unique(Location))) %>%
  mutate(MTV = median(StdValue), LogMTV = median(LogTraitValue), 
         Mean = mean(StdValue), LogMean = mean(LogTraitValue),
         SD = sd(StdValue), LogSD = sd(LogTraitValue), 
         COV = SD/Mean, LogCOV = LogSD/LogMean, 
         LogCOVnew = log(SD/Mean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, TraitShort, NumbRecords, NumbSites, MTV, LogMTV, SD, LogSD, COV, LogCOV, LogCOVnew)

# Height with 5 random records
hei.itv.random <- itv.db %>% filter(TraitShort == "PlantHeight") %>% 
  group_by(sp) %>% 
  slice_sample(n = 5) %>%
  mutate(SampleNumbRecords = 5, SampleNumbSites = length(unique(Location))) %>%
  mutate(SampleMTV = median(StdValue), LogSampleMTV = median(LogTraitValue), 
         SampleMean = mean(StdValue), SampleLogMean = mean(LogTraitValue),
         SampleSD = sd(StdValue), LogSampleSD = sd(LogTraitValue), 
         SampleCOV = SampleSD/SampleMean, SampleLogCOV = LogSampleSD/SampleLogMean,
         SampleLogCOVnew = log(SampleSD/SampleMean)) %>% 
  ungroup() %>% distinct(sp, .keep_all = TRUE) %>% 
  select(sp, SampleNumbRecords, SampleNumbSites, SampleMTV, LogSampleMTV, 
         SampleSD, LogSampleSD, SampleCOV, SampleLogCOV, SampleLogCOVnew)

# Full table
hei.itv.all <- left_join(hei.itv, hei.itv.random, by = "sp")

write.csv(hei.itv.all, "data/2022_height_itv_comparison.csv")




## ITV SCATTERPLOTS ----

# SLA scatterplot
(itv.sp.sla.newcov <- ggplot(sla.itv.all) + 
    geom_point(aes(x = LogSD, y = LogCOVnew), size = 4, colour = "#0a6933", alpha = 0.8) + 
    ylab("Log COV\n") + xlab("\nLog SD") +
    ylim(-2.5, 0.2) +
    theme(axis.text.x  = element_text(size = 16, colour = "black"), 
          legend.title = element_text(size = 25, face = "bold"), 
          legend.text=element_text(size = 25),
          axis.title.x = element_text(face="bold", size=22),
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=20, colour = "black"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          plot.margin = unit(c(3,2,2,2), "lines"),
          legend.background = element_blank(), legend.key = element_blank()))


# Seed mass scatterplot
(itv.sp.seed.newcov <- ggplot(seed.itv.all) + 
    geom_point(aes(x = LogSD, y = LogCOVnew), size = 4, colour = "#E57E00", alpha = 0.8) + 
    ylab("Log COV\n") + xlab("\nLog SD") +
    theme(axis.text.x  = element_text(size = 16, colour = "black"), 
          legend.title = element_text(size = 25, face = "bold"), 
          legend.text=element_text(size = 25),
          axis.title.x = element_text(face="bold", size=22),
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=20, colour = "black"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          plot.margin = unit(c(3,2,2,2), "lines"),
          legend.background = element_blank(), legend.key = element_blank()))


# Height scatterplot
(itv.sp.hei.newcov <- ggplot(hei.itv.all) + 
    geom_point(aes(x = LogSD, y = LogCOVnew), size = 4, colour = "#800080", alpha = 0.8) +
    ylab("Log COV\n") + xlab("\nLog SD") +
    theme(axis.text.x  = element_text(size = 16, colour = "black"), 
          legend.title = element_text(size = 25, face = "bold"), 
          legend.text=element_text(size = 25),
          axis.title.x = element_text(face="bold", size=22),
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=20, colour = "black"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          plot.margin = unit(c(3,2,2,2), "lines"),
          legend.background = element_blank(), legend.key = element_blank()))


## Panel (Figure S6)
(itv.sp.newcov.panel <- plot_grid(itv.sp.sla.newcov, itv.sp.seed.newcov, itv.sp.hei.newcov,  
                                  ncol=3, nrow = 1, align="hv", label_size = 26,
                                  labels = c("a) SLA", "b) Seed Mass", "c) Height")))

ggplot2::ggsave(itv.sp.newcov.panel, filename = "figures/Figure_S6.png", 
                width = 40, height = 20, units = "cm")

