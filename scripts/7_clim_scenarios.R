## Trait-range manuscript
## Mariana Garcia
## Jan 2021
## Script 7. Comparison of climate change scenarios


# Some input files for this script are not available in the repo given that they include raw data.
# The script preparing the data is included below for transparency and reproducibility.
# The full mastersheet (master_new_superfinal00.RData) is available below as input data for all analyses and figures in this script.


## LIBRARIES ----
library(tidyverse)
library(brms)
library(ggpubr)


## THEME ----
range.theme <- theme(legend.position = "none",
                     axis.title.x = element_text(face="bold", size=22),
                     axis.text.x  = element_text(vjust=0.5, size=20, colour = "black"), 
                     axis.title.y = element_text(face="bold", size=22),
                     axis.text.y  = element_text(vjust=0.5, size=20, colour = "black"),
                     panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
                     panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                     plot.title = element_text(color = "black", size = 18, face = "bold", hjust = 0.5),
                     plot.margin = unit(c(1,1,1,1), units = , "cm"))


## We need to re-make our mastersheet to incorporate all 3 climatic scenarios


## LOAD DATA ----

# Trait data with categories
load("data/2022_trait_cat.RData")
#this is assigned to the object 'trait.cat' - 17921 obs 

# Raw trait data
load("data/2022_try5_clean5.RData") # [this file is not available as it contains raw data]
name.vector.try5 <- unique(try5.clean5$sp)


## RANGE DATA ----

# [this file is not available in the repo as it contains raw data]
range.data2 <- read.csv("data/SRC_baseline_tabs_2017-04-26_fullrange.csv")
unique(range.data2$sp) #132 species

# Convert species names with dot to spaces to match TRY
range.data2$sp <- gsub(".", " ", range.data2$sp, fixed = TRUE)

# convert the two problematic hiphenated species
range.data2$sp[range.data2$sp == "Arctostaphylos uvaursi"] <- "Arctostaphylos uva-ursi"
range.data2$sp[range.data2$sp == "Vaccinium vitisidaea"] <- "Vaccinium vitis-idaea"

# Species name changes: convert species that have PAF nomenclature to The Plant List nomenclature
# I have done this manually comparing Word docs 
range.data2$sp[range.data2$sp == "Betula fruticosa"] <- "Betula humilis"
range.data2$sp[range.data2$sp == "Chamaepericlymenum canadense"] <- "Cornus canadensis"
range.data2$sp[range.data2$sp == "Swida sericea"] <- "Cornus sericea"
range.data2$sp[range.data2$sp == "Kalmia procumbens"] <- "Loiseleuria procumbens"
range.data2$sp[range.data2$sp == "Oxycoccus palustris"] <- "Vaccinium oxycoccos"


# Filter right parameters and convert pixel number into km (each pixel size is 10x10km)
# Divide range changes by 100 so the values are smaller and log-transforming works properly
range.all <- range.data2 %>% 
  filter(area == "full_area") %>% 
  filter(biointer == "no") %>%
  filter(sp %in% name.vector.try5) %>%
  mutate(scenario = paste0(rcp, "_", gcm)) %>%
  mutate(RangeLog = log(CurrentRangeSize)) %>%
  mutate(CurrentRangeSizeKm = CurrentRangeSize*100) %>%
  mutate(RangeLogKm = log(CurrentRangeSizeKm)) %>%
  mutate(FutureRangeSize.FullDisp.Km = FutureRangeSize.FullDisp*100) %>% 
  mutate(AbsoluteRangeChangeKm = FutureRangeSize.FullDisp.Km - CurrentRangeSizeKm) %>%
  mutate(AbsoluteRangeChangeKmDivided = AbsoluteRangeChangeKm/1000000) %>%
  mutate(SpeciesRangeChangeDivided = SpeciesRangeChange/100)

# Find the smallest value so we can use it to add as a constant
min.abs <- range.all %>% dplyr::select(AbsoluteRangeChangeKmDivided) %>% min() %>% abs()
min.rel <- range.all %>% dplyr::select(SpeciesRangeChangeDivided) %>% min() %>% abs()
min.test <- range.all %>% dplyr::select(SpeciesRangeChange) %>% min() %>% abs()
min.abs.test <- range.all %>% dplyr::select(AbsoluteRangeChangeKm) %>% min() %>% abs()

# Add constant so we can log-transform negative values, calculate quantiles, log-transform and centre data
range.all2 <- range.all %>%
  dplyr::mutate(AbsoluteRangeChangeKmConstant = AbsoluteRangeChangeKmDivided + min.abs + 1) %>%
  dplyr::mutate(SpeciesRangeChangeConstant = SpeciesRangeChangeDivided + min.rel + 1) %>%
  dplyr::group_by(sp, filt) %>% 
  dplyr::mutate(abs.median = quantile(AbsoluteRangeChangeKmConstant, probs = 0.50)) %>%
  ungroup() %>%
  dplyr::mutate(abs.median.log = log(abs.median)) %>%
  dplyr::distinct(sp, filt, .keep_all = TRUE) %>%
  dplyr::group_by(filt) %>% 
  dplyr::mutate(cent.abs.median.log = abs.median.log - mean(abs.median.log))

# merge with traits
merge.df <- merge(try5.clean5, range.all2, by = "sp")


## Add decidiousness 

# [this file is not available in the repo as it contains raw data]
form <- read.csv("data/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.csv")
vector.xx <- unique(merge.df$sp)

# filter for the relevant species in our list
form.data <- form %>% filter(AccSpeciesName %in% vector.xx)
form.short <- dplyr::select(form.data, AccSpeciesName, Family, LeafPhenology)
form.short$AccSpeciesName # data for 54 species

# join with main trait database
names(form.short)[names(form.short) == "AccSpeciesName"] <- "sp"
master.full <- left_join(merge.df, form.short, by = "sp")

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

# Modify with capitals for the graphs
master.full$LeafPhenology <- as.character(master.full$LeafPhenology)
master.full$LeafPhenology[master.full$LeafPhenology == "evergreen"] <- "Evergreen"
master.full$LeafPhenology[master.full$LeafPhenology == "deciduous"] <- "Deciduous"


## Fix functional groups
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


## Dispersal

# Collect berry species in a vector
berry.vector <- c("Vaccinium caespitosum", "Vaccinium myrtilloides", "Vaccinium myrtillus", 
                  "Vaccinium uliginosum", "Vaccinium vitis-idaea", "Vaccinium oxycoccos", 
                  "Empetrum nigrum", "Rubus chamaemorus",
                  "Rubus idaeus", "Cornus canadensis", "Arctostaphylos uva-ursi", "Arctous alpina",
                  "Arctous rubra", "Ribes spicatum", "Ribes glandulosum", "Ribes nigrum")

# Add in dispersal mode per species
master.new.superfinal00 <- mutate(master.full, DispersalMode = ifelse(sp %in% berry.vector, "Berry", "Wind"))

# Keep relevant columns only
master.new.superfinal00 <- master.new.superfinal00 %>% 
  select(-c(X, layer_id, Loss, Stable0, Stable1, Gain, PercLoss, PercGain, SpeciesRangeChange, 
            CurrentRangeSize, CurrentRangeSizeKm, FutureRangeSize.FullDisp, FutureRangeSize.FullDisp, FutureRangeSize.FullDisp.Km,
            area, file.id, model, rcp, gcm, src_ras_file, scenario, RangeLog, RangeLogKm, FutureRangeSize.FullDisp.Km, 
            AbsoluteRangeChangeKm, AbsoluteRangeChangeKmDivided, AbsoluteRangeChangeKmConstant, SpeciesRangeChangeDivided,
            SpeciesRangeChangeConstant, abs.median, abs.median.log))

save(master.new.superfinal00, file = "data/master_new_superfinal00.RData")




## SCENARIOS ----

# This file can be loaded to start fitting the models below
load("data/master_new_superfinal00.RData")

unlim <- filter(master.new.superfinal00, filt == "unlimited_dipersal")
lim.disp <- filter(master.new.superfinal00, filt == "max_dipersal")
no.disp <- filter(master.new.superfinal00, filt == "no_dipersal")

## Limited dispersal - seed
limdisp.seed.sp <- lim.disp %>%  
  filter(TraitShort == "SeedMass") %>% 
  mutate(LogTraitValue = log(StdValue)) %>% 
  filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  group_by(sp) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE)

limdisp.seed.sp.mod <- brm(cent.abs.median.log ~ LogTraitMedianValue, data = limdisp.seed.sp, 
                          iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_lim_disp_seed_mod")
summary(limdisp.seed.sp.mod) # negative ns
bayes_R2(limdisp.seed.sp.mod) #R2 = 4.9


## Unlimited dispersal - seed
unlimdisp.seed.sp <- unlim %>%  
  filter(TraitShort == "SeedMass") %>% 
  mutate(LogTraitValue = log(StdValue)) %>% 
  filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  group_by(sp) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE)

unlim.seed.sp.mod <- brm(cent.abs.median.log ~ LogTraitMedianValue, data = unlimdisp.seed.sp, 
                           iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_unlim_disp_seed_mod")
summary(unlim.seed.sp.mod) # negative ns
bayes_R2(unlim.seed.sp.mod) #R2 = 4.3




## Limited dispersal - height
limdisp.hei.sp <- lim.disp %>%  
  filter(TraitShort == "PlantHeight") %>% 
  mutate(LogTraitValue = log(StdValue)) %>% 
  filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  group_by(sp) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE)

limdisp.hei.sp.mod <- brm(cent.abs.median.log ~ LogTraitMedianValue, data = limdisp.hei.sp, 
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_lim_disp_hei_mod")
summary(limdisp.hei.sp.mod) # positive ns
bayes_R2(limdisp.hei.sp.mod) #R2 = 4.4%


## Unlimited dispersal - height
unlim.hei.sp <- unlim %>%  
  filter(TraitShort == "PlantHeight") %>% 
  mutate(LogTraitValue = log(StdValue)) %>% 
  filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  group_by(sp) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE)

unlim.hei.sp.mod <- brm(cent.abs.median.log ~ LogTraitMedianValue, data = unlim.hei.sp, 
                          iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_unlim_disp_hei_mod")
summary(unlim.hei.sp.mod) # positive ns
bayes_R2(unlim.hei.sp.mod) #R2 = 2.37%



## Unlimited dispersal (climatic scenario only)

# Height data
unlim.hei.sp.wei <- unlim %>%  
  dplyr::filter(TraitShort == "PlantHeight") %>% 
  dplyr::mutate(LogTraitValue = log(StdValue)) %>% 
  dplyr::filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  dplyr::group_by(sp) %>% 
  dplyr::mutate(TraitValueSD = sd(StdValue)) %>%
  dplyr::mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  dplyr::mutate(MedianTraitValue = median(StdValue)) %>%
  dplyr::mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>% 
  dplyr::mutate(index = ifelse(nobs >=20, 1,
                        ifelse(nobs < 20, 0.33+(nobs/30), 0)))

# Model
unlim.hei.sp.wei.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                            data = unlim.hei.sp.wei, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/2022_unlim_disp_wei_hei_mod")
summary(unlim.hei.sp.wei.mod) # positive ns



# Predictions
heivalabs20 <- data.frame(LogTraitMedianValue = 
                          seq(from = min(unlim.hei.sp.wei$LogTraitMedianValue),
                              to = max(unlim.hei.sp.wei$LogTraitMedianValue), by = 0.1))
# predict values
heifit600 <- fitted(
  unlim.hei.sp.wei.mod, 
  newdata = heivalabs20, 
  re_formula = NULL,
  summary = TRUE
)

# combine dataframes
colnames(heifit600) = c('fit', 'se', 'lwr', 'upr')
df_hei600 = cbind(heivalabs20, heifit600)


## Plot relationships
(hei.val.medabs.unlim.plot <- ggplot() + 
    geom_point(data = unlim.hei.sp.wei, aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                                   colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nHeight values (log m)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_line(data = df_hei600, aes(x = LogTraitMedianValue, y = fit), colour = "#800080") + 
    geom_ribbon(data = df_hei600, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), fill = "#800080", alpha = 0.3) +
    range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE)))



#### SLA data
unlim.sla.sp.wei <- unlim %>%  
  dplyr::filter(TraitShort == "SLA") %>% 
  dplyr::mutate(LogTraitValue = log(StdValue)) %>% 
  dplyr::filter_at(vars(LogTraitValue), all_vars(!is.infinite(.))) %>%
  dplyr::group_by(sp) %>% 
  dplyr::mutate(TraitValueSD = sd(StdValue)) %>%
  dplyr::mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  dplyr::mutate(MedianTraitValue = median(StdValue)) %>%
  dplyr::mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>% 
  dplyr::mutate(index = ifelse(nobs >=20, 1,
                               ifelse(nobs < 20, 0.33+(nobs/30), 0)))

# Model
unlim.sla.sp.wei.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = unlim.sla.sp.wei, 
                        iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_unlim_disp_wei_sla_mod")
summary(unlim.sla.sp.wei.mod) #negative ns


# Predictions
slavalabs20 <- data.frame(LogTraitMedianValue = 
                            seq(from = min(unlim.sla.sp.wei$LogTraitMedianValue),
                                to = max(unlim.sla.sp.wei$LogTraitMedianValue), by = 0.1))
# predict values
slafit600 <- fitted(
  unlim.sla.sp.wei.mod, 
  newdata = slavalabs20, 
  re_formula = NULL,
  summary = TRUE
)

# combine dataframes
colnames(slafit600) = c('fit', 'se', 'lwr', 'upr')
df_sla600 = cbind(slavalabs20, slafit600)



## Plot relationships
(sla.val.medabs.unlim.plot <- ggplot() + 
    geom_point(data = unlim.sla.sp.wei, aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                                        colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("\nSLA values (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Deciduousness") +
    geom_line(data = df_sla600, aes(x = LogTraitMedianValue, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_sla600, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), fill = "#0a6933", alpha = 0.5) +
    range.theme + 
    theme(legend.title = element_text(size = 20, face = "bold"), 
        legend.text=element_text(size = 20),
        legend.position = "top", legend.key = element_blank(),
        legend.background = element_blank()))




## Seed Mass
unlim.seed.sp <- unlim %>% 
  dplyr::filter(TraitShort == "SeedMass") %>% 
  dplyr::mutate(StdValue = replace(StdValue, StdValue == 0, 0.0001)) %>%
  dplyr::mutate(LogTraitValue = log(StdValue)) %>%
  dplyr::group_by(sp) %>% 
  dplyr::mutate(TraitValueSD = sd(StdValue)) %>%
  dplyr::mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  dplyr::mutate(MedianTraitValue = median(StdValue)) %>%
  dplyr::mutate(LogTraitMedianValue = median(LogTraitValue)) %>% 
  dplyr::distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>% 
  dplyr::mutate(index = ifelse(nobs >=20, 1,
                        ifelse(nobs < 20, 0.33 +(nobs/30), 0)))



## Additional step in here as we need to bring in the gap-filled data
gap.seeds.fut <- read.csv("data/2022_gapfilled_seed_data.csv")

# extract unique columns so they can be combined
gap.unique <- gap.seeds.fut %>% dplyr::select(sp, GapMedian, GapSD)

# extract the species vector with the gap-filled species
gap.vector <- unique(gap.seeds.fut$sp)

# extract all the range columns from the principal database above for the gap-filled species
gap.long <- unlim %>% filter(sp %in% gap.vector) %>% distinct(sp, .keep_all = TRUE)

# merge into one dataframe
gap.more <- merge(gap.long, gap.unique, by = "sp")

# not sure why the distinct() is not working above - shoudl be 12 obs
gap.more.more <- distinct(gap.more, sp, .keep_all = TRUE)

# gap.seeds.fut needs to have the same columns as seed.quan so they can be binded
# we add an index value of 0.5 as we trust less the species that we have been gap-filling
# remove redundant columns to facilitate row bindng
# I am replacing the 0 values of the SD by half of the minimum value that is not zero
gap.seeds.fut.new <- gap.more.more %>% 
  mutate(LogTraitValue = NA) %>% 
  mutate(TraitValueSD = NA) %>%
  mutate(LogTraitValueSD = GapSD) %>% 
  mutate(MedianTraitValue = NA) %>%
  mutate(LogTraitMedianValue = GapMedian) %>% 
  mutate(index = 0.5) %>%
  mutate(TraitShort = "SeedMass") %>%
  mutate(Dataset = "Gap-filled") %>% 
  mutate(ValueKindName = "Gap-filled") %>% mutate(DataContributor = "Gap-filled") %>% 
  mutate(StdValue = NA) %>% mutate(Source = "Gap-filled") %>%
  dplyr::select(., -c(GapMedian, GapSD))


# reorder columns in gap.seeds.new according to column order in seed.now
gap.seeds.fut.new2 <- gap.seeds.fut.new[names(unlim.seed.sp)]

# bind dataframes into one
unlim.seed.all <- bind_rows(unlim.seed.sp, gap.seeds.fut.new2)
str(unlim.seed.all)



# Model
unlim.seed.sp.wei.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = unlim.seed.all, 
                        iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_unlim_disp_wei_seed_mod")
summary(unlim.seed.sp.wei.mod) #negative ns



# Predictions
seedvalabs20 <- data.frame(LogTraitMedianValue = 
                            seq(from = min(unlim.seed.all$LogTraitMedianValue),
                                to = max(unlim.seed.all$LogTraitMedianValue), by = 0.1))
# predict values
seedfit600 <- fitted(
  unlim.seed.sp.wei.mod, 
  newdata = seedvalabs20, 
  re_formula = NULL,
  summary = TRUE
)

# combine dataframes
colnames(seedfit600) = c('fit', 'se', 'lwr', 'upr')
df_seed600 = cbind(seedvalabs20, seedfit600)


# Let's filter depending on source
gap.points.fut.x <- filter(unlim.seed.all, Source == "Gap-filled")
proper.sources.x <- c("LedaKleyer", "TTT", "TRY", "BobHollisterSeeds", "LucieSmrzova", "EstherLevesqueBylot", "No")
other.points.fut.x <- filter(unlim.seed.all, Source %in% proper.sources.x)

# Plotting relationships
(seed.val.unlim.plot <- ggplot() + 
    geom_point(data = other.points.fut.x,
               aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points.fut.x, shape = 1, 
               aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    xlab("\nSeed Mass values (log mg)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seed600, aes(x = LogTraitMedianValue, y = fit), colour = "#E57E00") + 
    geom_ribbon(data = df_seed600, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr),
                fill = "#E57E00", alpha = 0.4) +
    labs(colour = "Dispersal mode") +
    range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))



## Height variation
hei.var.unlim.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = unlim.hei.sp.wei,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_unlim_disp_wei_hei_var_mod")
summary(hei.var.unlim.mod) # positive ns


# Model predictions
hei.var.unlim.data = data.frame(LogTraitValueSD = seq(from = min(unlim.hei.sp.wei$LogTraitValueSD),
                                                to = max(unlim.hei.sp.wei$LogTraitValueSD), by = 0.1))

fit.heivarunlim = fitted(
  hei.var.unlim.mod,
  newdata = hei.var.unlim.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.heivarunlim) = c('fit', 'se', 'lwr', 'upr')
df_heivarunlim = cbind(hei.var.unlim.data, fit.heivarunlim)

# Plot relationships
(hei.var.unlim.plot <- ggplot() + 
    geom_point(data = unlim.hei.sp.wei, aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                                   colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nHeight variation (log m)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_line(data = df_heivarunlim, aes(x = LogTraitValueSD, y = fit), colour = "#800080") + 
    geom_ribbon(data = df_heivarunlim, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), fill = "#800080", alpha = 0.3) +
    range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE)))



## SLA variation
sla.unlimvar.wei.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = unlim.sla.sp.wei,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_unlim_disp_wei_sla_var_mod")
summary(sla.unlimvar.wei.mod) # positive ns



## Model predictions
sla.unlimvar.data = data.frame(LogTraitValueSD = seq(from = min(unlim.sla.sp.wei$LogTraitValueSD),
                                                to = max(unlim.sla.sp.wei$LogTraitValueSD), by = 0.1))

fit.unlimvar.sla = fitted(
  sla.unlimvar.wei.mod,
  newdata = sla.unlimvar.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.unlimvar.sla) = c('fit', 'se', 'lwr', 'upr')
df_sla_unlimvar = cbind(sla.unlimvar.data, fit.unlimvar.sla)


## Plot relationships
(sla.unlimvar.plot <- ggplot() + 
    geom_point(data = unlim.sla.sp.wei, aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                                   colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("\nSLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_sla_unlimvar, aes(x = LogTraitValueSD, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_sla_unlimvar, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), 
                fill = "#0a6933", alpha = 0.5) + range.theme + 
    labs(colour = "Deciduousness") +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))


## Seed Mass variation
seed.unlimvar.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD,
                       data = unlim.seed.all, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_unlim_disp_wei_seed_var_mod")
summary(seed.unlimvar.mod) # positive ns



# Model predictions
seed.unlimvar.data = data.frame(LogTraitValueSD = seq(from = min(unlim.seed.all$LogTraitValueSD),
                                                 to = max(unlim.seed.all$LogTraitValueSD), by = 0.1))

fit.seed.unlimvar = fitted(
  seed.unlimvar.mod,
  newdata = seed.unlimvar.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.seed.unlimvar) = c('fit', 'se', 'lwr', 'upr')
df_seed_unlimvar = cbind(seed.unlimvar.data, fit.seed.unlimvar)


# Plot relationships
(seed.unlimvar.plot <- ggplot() + 
    geom_point(data = other.points.fut.x,
               aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points.fut.x, shape = 1, 
               aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    xlab("\nSeed Mass Variation (log mg)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seed_unlimvar, aes(x = LogTraitValueSD, y = fit), colour = "#E57E00") + 
    geom_ribbon(data = df_seed_unlimvar, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr),
                fill = "#E57E00", alpha = 0.4) +
    labs(colour = "Dispersal mode") +
    range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))


## PANEL
(unlim.panel <- ggarrange(hei.val.medabs.unlim.plot, sla.val.medabs.unlim.plot, seed.val.unlim.plot, 
                          hei.var.unlim.plot, sla.unlimvar.plot, seed.unlimvar.plot,
                            labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
                            nrow = 2, ncol = 3, font.label = list(size = 30)))
ggplot2::ggsave(unlim.panel, filename = "figures/Figure_S2.png", 
                width = 60, height = 60, dpi = 500, units = "cm")




## COMPARISON OF CLIMATIC SCENARIOS ----
clim.scen <- range.data2 # [this file is not available in the repo as it contains raw range data]

# Filter for appropriate parameters and calculate median range change per climatic scenario and species
clim.scen.fin <- clim.scen %>% filter(area == "full_area") %>% filter(biointer == "no") %>% 
  filter(sp %in% name.vector.try5) %>% mutate(CurrentRangeSizeKm = CurrentRangeSize*100) %>%
  mutate(FutureRangeSize.FullDisp.Km = FutureRangeSize.FullDisp*100) %>% 
  mutate(AbsoluteRangeChangeKm = FutureRangeSize.FullDisp.Km - CurrentRangeSizeKm) %>%
  group_by(sp, filt) %>% mutate(MedianRangeChange = median(AbsoluteRangeChangeKm)) %>%
  distinct(sp, filt, .keep_all = TRUE)

# Order in descending range change for the limited dispersal scenario
limited.order <- clim.scen.fin %>% filter(filt == "max_dipersal") %>% arrange(-MedianRangeChange)
limited.order.vector <- unique(limited.order$sp)

# Plot the values per species and climatic scenario (Fig S1)
(clim.scen.plot <- ggplot(clim.scen.fin, aes(x = factor(sp, level = limited.order.vector), y = MedianRangeChange, colour = filt)) + 
    geom_point(size = 6) + 
    geom_hline(yintercept = 0, linetype = "solid") +
    xlab("\nSpecies") + ylab(expression(bold(paste("Absolute Species Range Change (km"^bold("2"), ")")))) +
    scale_colour_manual(values = c("#44BB7E", "#EDBD12", "#9A39C6"), labels=c("Limited Dispersal", "No Dispersal", "Unlimited Dispersal")) + 
    scale_y_continuous(breaks = seq(-7850000, 27288000, 5000000),
                       labels = function(x) format(x, scientific = FALSE)) +
    labs(colour = "Climatic scenario") + range.theme + 
    theme(axis.text.x  = element_text(face = "italic", angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          axis.title.y = element_text(face="bold", size=22),
          legend.title = element_text(size = 22, face = "bold"), 
          legend.text=element_text(size = 22),
          legend.position = c(0.9, 0.7), legend.key = element_blank(),
          legend.background = element_blank()))

ggsave(clim.scen.plot, filename = "figures/Figure_S1.png", width = 50, height = 30, units = "cm")

