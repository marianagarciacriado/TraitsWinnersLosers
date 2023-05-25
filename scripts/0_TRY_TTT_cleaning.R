### CLEAN TRY SPECIES NAMES ###
# Match to The Plant List
# Anne Bjorkman with further edits from Mariana
# March 2020


# The input files for this script are not available in the repo given that they include raw TRY/TTT data.
# This is a data prep script and is made available for reproducibility and transparency purposes. 
# The final file produced by this script (mTTTT_clean.RData) will be used as input data on species traits 
# in script #1 in order to produce the main mastersheet combining traits and ranges.



## PACKAGES ----
library(data.table)
library(dplyr)
library(Taxonstand)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)



## IMPORT DATA ----

# [this file is not available in the repo as it contains raw data]
TRYv5 <- fread("data/8400.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T, verbose=TRUE)
head(TRYv5)

allTRYnames <- unique(TRYv5$AccSpeciesName)


# Add location info [this file is not available in the repo as it contains raw data]
siteinfo <- read.csv(file = "data/TRY_5_Site_Climate_Soil_2019-03-25.txt", 
                     header = T, sep = "\t", dec = ".", stringsAsFactors=FALSE, strip.white=TRUE)

TRYv5$Lat <- siteinfo$LAT_site[match(TRYv5$ObservationID,siteinfo$observationId)]
TRYv5$Lon <- siteinfo$LON_site[match(TRYv5$ObservationID,siteinfo$observationId)]

# Import Mariana's list of species names
# SEE FILE: Mariana_TRY_species_list.txt

mnames <- c(name.vector, syn.vector, try50.vector)

mTRY <- as.data.frame(TRYv5[TRYv5$AccSpeciesName %in% c(mnames),], stringsAsFactors == F)


# SAVE this so you can just load this file and not the whole 19 GB TRY database if you want!
# save(mTRY, file = "data/TRYv5_Mariana_unclean.RData")



# START HERE IF JUST USING THE DATA SUBSET! ----

# [this file is not available in the repo as it contains raw data]
load("data/TRYv5_Mariana_unclean.RData") #object is called "mTRY"

## Checks on duplicates here

# Alastair Rogers
alst <- mTRY %>% filter(LastName == "Rogers") %>%
  filter(TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded")

unique(alst$TraitName)

# there are duplicates here
write.csv(alst, file = "data/alastair_rogers_b4_try.csv") # [this file is not available in the repo as it contains raw data]


# Serge Sheremetev
serge <- mTRY %>% filter(LastName == "Sheremetev") %>% filter(OriglName == "SLA") %>%
  filter(AccSpeciesName == "Vaccinium myrtillus")

# same here
write.csv(serge, file = "data/serge_sheremetev_b4_try.csv") # [this file is not available in the repo as it contains raw data]



# CHECK NAMES IN THE PLANT LIST ----
TRYcorr <- TPL(unique(mTRY$AccSpeciesName))
head(TRYcorr)

TRYcorr[TRYcorr$New.Species != TRYcorr$Species,]
TRYcorr[TRYcorr$Typo == "TRUE",]
unique(TRYcorr$New.Taxonomic.status)
TRYcorr[TRYcorr$Plant.Name.Index==FALSE,]

# All names seem to be accepted PL names (so they should match with TTT below, which has already been checked with the PL)

# Add Family & Genus data to dataframe (for cleaning)
mTRY$Genus <- TRYcorr$Genus[match(mTRY$AccSpeciesName,TRYcorr$Taxon)]
mTRY$Family <- TRYcorr$Family[match(mTRY$AccSpeciesName,TRYcorr$Taxon)]



# SUBSET TRY DATA (data type) ----

unique(mTRY$TraitName)
unique(mTRY$ValueKindName)
unique(mTRY$UncertaintyName)
unique(mTRY$DataID)

# Keep only the traits Mariana needs: Seed mass, plant height and SLA (all versions)
unique(mTRY$TraitName)[grep('SLA', unique(mTRY$TraitName))]
unique(mTRY$TraitName)[grep('height', unique(mTRY$TraitName))]
unique(mTRY$TraitName)[grep('Seed', unique(mTRY$TraitName))]

mTRY <- mTRY[mTRY$TraitName %in% c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded","Leaf area per leaf fresh mass (specific leaf area (SLA or 1/LMA) based on leaf fresh mass)","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded","Plant height vegetative","Plant height generative","Seed dry mass"),]

# Keep only data that are NOT from an "exposition" (DataID 327)

mTRY[mTRY$DataID == 327,] #none



# ADD TTT DATA ----

load("data/try_ttt.RData") # [this file is not available in the repo as it contains raw data]

# Let's filter species that I want only from TTT
try.ttt2 <- as.data.frame(try.ttt[try.ttt$AccSpeciesName %in% c(mnames),])
unique(try.ttt2$AccSpeciesName)

# Prepare TRY data to match with TTT
mTRY$DataContributor <- paste(mTRY$FirstName,mTRY$LastName, sep=" ")

mTRY$TraitShort <- mTRY$TraitName
mTRY$TraitShort[mTRY$TraitShort %in% c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "Leaf area per leaf fresh mass (specific leaf area (SLA or 1/LMA) based on leaf fresh mass)","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included","Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded")] <- "SLA"
mTRY$TraitShort[mTRY$TraitShort %in% c("Plant height vegetative","Plant height generative")] <- "PlantHeight"
mTRY$TraitShort[mTRY$TraitShort == "Seed dry mass"] <- "SeedMass"

mTRY$Source <- "TRY"
mTRY$Treatment <- NA

try.ttt2$ObservationID <- paste("TTT", seq(from = 1, to = nrow(try.ttt2)), sep="_") 

try.ttt2$ObsDataID <- paste("TTT", seq(from = 1, to = nrow(try.ttt2)), sep="_") # NOTE that this will be a unique number even for traits measured on the same individual. If you need to know which traits were measured on the same individual, that information is available in the TraitHub (data paper) data version of TTT

mTTTT <- rbind(mTRY[,intersect(colnames(mTRY), colnames(try.ttt2))], 
               try.ttt2[try.ttt2$Source == "TTT" & try.ttt2$TraitShort %in% c("PlantHeight","SeedMass","SLA"),])


# let's check duplicates here on TTT
ttt.test <- mTTTT %>% filter(Source == "TTT") %>% 
  group_by(AccSpeciesName, TraitShort, Lat, Lon, StdValue) %>% 
  mutate(dupe = n()>1) %>% filter(dupe == "TRUE") # 12886 obs

ttt.test.summ <- ttt.test %>% 
  group_by(TraitShort, Dataset, DataContributor) %>% 
  mutate(nobs = length(StdValue))

summ.ttt <- count(ttt.test.summ, Dataset)



# CLEAN COMBINED DATA ----

nrow(mTTTT)

# REMOVE NIWOT MEAN DATA (because individual values are added from TTT) and ALL NIWOT SLA DATA because of inconsistency
mTTTT <- mTTTT[c(mTTTT$Dataset == "Niwot Alpine Plant Traits" | c(mTTTT$Dataset == "NiwotRidge" & mTTTT$TraitShort == "SLA"))==F, ]

# REMOVE SEED MASS DATA from Salix arctica Rebecca Klady because units uncertain and Papaver from me because values seem off (units problem?)
mTTTT <- mTTTT[c(mTTTT$DataContributor=="Anne Bjorkman" & mTTTT$AccSpeciesName=="Papaver radicatum" & mTTTT$TraitShort=="SeedMass")==F,]
mTTTT <- mTTTT[c(mTTTT$DataContributor=="Rebecca Klady" & mTTTT$AccSpeciesName=="Salix arctica" & mTTTT$TraitShort=="SeedMass")==F,]


# REMOVE ALL TREATMENT DATA AND BOTANICAL GARDENS

unique(mTTTT$Treatment)
mTTTT <- mTTTT[mTTTT$Treatment %in% c("none","control","Natural environment","None","Control","Warming over time") | is.na(mTTTT$Treatment),]

# KEEP ONLY ESTIMATES THAT ARE MEASURED ON INDIVIDUALS OR APPROXIMATE THE MEAN
unique(mTTTT$ValueKindName)
mTTTT <- mTTTT[mTTTT$ValueKindName %in% c("Maximum","Minimum","Low","High","Maximum in plot") == F,]



# AUTOMATED CLEANING ----

# Add genus name
mTTTT$Genus<-as.vector(sapply(strsplit(mTTTT$AccSpeciesName," ",fixed=FALSE), "[", 1))

(nrow.start <- nrow(mTTTT))
# 20218 obs so far

ttt.unclean <- mTTTT

## Step 1: individual records against the entire distribution (except plant height and leaf area)
ttt.clean2 <- ttt.unclean %>% dplyr::mutate(nrow=row_number()) %>% group_by(TraitShort) %>% dplyr::mutate(
  mean_all = round(((sum(StdValue) - StdValue)/(n()-1)),5), 
  n_all = n(), 
  median_all = median(StdValue),
  sd_all = sd(StdValue), 
  ErrorRisk_all = round((abs(StdValue-mean_all)/sd_all),4),
  ErrorRiskMedian_all = round((abs(StdValue-median_all)/sd_all),4))

ttt.clean3 <- ttt.clean2[(is.finite(ttt.clean2$ErrorRisk_all) & ttt.clean2$ErrorRisk_all > 8 & ttt.clean2$TraitShort %in% c("PlantHeight") == F) == FALSE,]

## Step 2: datasets-by-species against species distribution
ttt.clean4 <- ttt.clean3 %>% group_by(TraitShort, AccSpeciesName) %>% dplyr::mutate(
  ndataset = length(unique(Dataset)), 
  mean_species = round(((sum(StdValue) - StdValue)/(n()-1)),5), 
  median_species = median(StdValue),
  sd_species = sd(StdValue), 
  n_species = n(), 
  ErrorRisk_species = round((abs(StdValue-mean_species)/sd_species),4),
  ErrorRiskMedian_species = round((abs(StdValue-median_species)/sd_species),4)) %>% 
  
  group_by(TraitShort, Genus) %>% dplyr::mutate(
    mean_genus = round(((sum(StdValue) - StdValue)/(n()-1)),5),
    median_genus = median(StdValue),
    sd_genus = sd(StdValue), 
    n_genus = n(), 
    ErrorRisk_genus = round((abs(StdValue-mean_genus)/sd_genus),4),
    ErrorRiskMedian_genus = round((abs(StdValue-median_genus)/sd_genus),4)) %>% 
  
  group_by(TraitShort, AccSpeciesName, Dataset) %>% dplyr::mutate(
    mean_spdataset = round(((sum(mean_species) - mean_species)/(n()-1)),5), 
    median_spdataset = median(StdValue),
    sd_spdataset = sd(StdValue), 
    n_spdataset = n(), 
    ErrorRisk_spdataset = round((abs(mean_spdataset-mean_species)/sd_species),4),
    ErrorRiskMedian_spdataset = round((abs(median_spdataset-median_species)/sd_species),4)) %>% 
  
  group_by(TraitShort, Genus, Dataset) %>% dplyr::mutate(
    mean_gspdataset = round(((sum(mean_genus) - mean_genus)/(n()-1)),5), 
    median_gspdataset = median(StdValue),
    sd_gspdataset = sd(StdValue), 
    n_gdataset = n(), 
    ErrorRisk_gspdataset = round((abs(mean_spdataset-mean_genus)/sd_genus),4),
    ErrorRiskMedian_gspdataset = round((abs(median_spdataset-median_genus)/sd_genus),4))

# > 4 datasets per species
ttt.clean5 <- ttt.clean4[(ttt.clean4$ndataset >= 4 & is.finite(ttt.clean4$ErrorRisk_spdataset) & ttt.clean4$ErrorRisk_spdataset > 3) == FALSE,] 

## Step 3: datasets-by-species against genus distribution
ttt.clean6 <- ttt.clean5[(ttt.clean5$ndataset >= 1 & ttt.clean5$ndataset < 4 & is.finite(ttt.clean5$ErrorRisk_gspdataset) & ttt.clean5$ErrorRisk_gspdataset > 3.5) == FALSE,]

## Step 4: individual records against the species distribution
ttt.clean7 <- ttt.clean6 %>% group_by(TraitShort, AccSpeciesName) %>% dplyr::mutate(
  ndataset = length(unique(Dataset)), 
  mean_species = round(((sum(StdValue) - StdValue)/(n()-1)),5), 
  mad_species = mad(StdValue, constant=1), 
  sd_species = sd(StdValue), n_species = n(), 
  ErrorRisk_species = round((abs(StdValue-mean_species)/sd_species),4))

# 1 - 3 records - leave everything as-is

# 4 - 9 records
ttt.clean8 <- ttt.clean7[(ttt.clean7$n_species >= 4 & ttt.clean7$n_species < 10 & is.finite(ttt.clean7$ErrorRisk_species) & ttt.clean7$ErrorRisk_species > 2.25) == FALSE,]

# 10 - 19 records
ttt.clean9 <- ttt.clean8[(ttt.clean8$n_species >= 10 & ttt.clean8$n_species < 20 & is.finite(ttt.clean8$ErrorRisk_species) & ttt.clean8$ErrorRisk_species > 2.75) == FALSE,]

# 20 - 29 records
ttt.clean10 <- ttt.clean9[(ttt.clean9$n_species >= 20 & ttt.clean9$n_species < 30 & is.finite(ttt.clean9$ErrorRisk_species) & ttt.clean9$ErrorRisk_species > 3.25) == FALSE,]

# >30 records
ttt.clean.final <- ttt.clean10[(ttt.clean10$n_species >= 30 & is.finite(ttt.clean10$ErrorRisk_species) & ttt.clean10$ErrorRisk_species > 4) == FALSE,]

mTTTT.clean <- as.data.frame(ttt.clean.final)

nrow.start - nrow(mTTTT.clean) #number of observations removed (186)


# SAVE CLEANED DATA
mTTTT.clean <- mTTTT.clean[,1:13]

save(mTTTT.clean, file = "data/mTTTT_clean.RData") # [this file is not available in the repo as it contains raw data]

unique(mTTTT.clean$AccSpeciesName)



# VIEW CLEANED DATA ----

# Check data removal
mTTTT.clean$NotRemoved <- 1
unmatched <- merge(ttt.unclean, mTTTT.clean, all=T)
unmatched$NotRemoved[is.na(unmatched$NotRemoved)] <- 0
unmatched$NotRemoved <- factor(unmatched$NotRemoved, levels = c(1,0))

# Histograms of kept vs. removed values per species and trait
# NOTE: this can take a long time to run

d_ply(unmatched, .(TraitShort,AccSpeciesName), function (x){
  ggplot(data=x)+
    geom_histogram(aes(StdValue,fill=NotRemoved),bins=200)+
    theme_bw()+
    ggtitle(paste(unique(x$AccSpeciesName),unique(x$TraitShort),sep="_") )
}, .print = T)

dev.off()

