## Trait-range manuscript
## Mariana Garcia
## August 2019
## Script 3. Current trait-range relationships


## LIBRARIES ----
library(brms)
library(tidyverse)
library(ggpubr)
library(modelr)
library(rstan)
library(ggrepel)



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


## LOAD & PREP DATA ----
load("data/2022_master_new_superfinal.RData")
#this is assigned to the object 'master.new.superfinal' - 17921 obs

# load winner/loser categories from script#4
cats <- read.csv("data/range_quan.csv")

# select only category
cats.short <- cats %>% dplyr::select(sp, category.abs)

# add category info on the mastersheet
trait.cat0 <- left_join(master.new.superfinal, cats.short, by = "sp")

# Add shortened species names
trait.cat0$SpCopy = trait.cat0$sp
trait.cat0 <- trait.cat0 %>% separate(SpCopy, c("GenShort", "SpShort"), " ") 
trait.cat <- trait.cat0 %>% mutate(GenShortName = str_sub(GenShort, start = 1, end = 3), 
                                  SpShortName = str_sub(SpShort, start = 1, end = 3)) %>% 
  mutate(SpLabel = paste0(GenShortName, " ", SpShortName))

# save for script#4
save(trait.cat, file = "data/2022_trait_cat.RData")

# Each row has an individual trait value, and we need we need only a line per species 
# (e.g. a single value of species' trait variation, median AND current range)




## Create objects by filtering per trait

# filter for SLA and log-transform the individual trait records
# group per species, calculate SD and median
sla.now <- trait.cat %>% 
  filter(TraitShort == "SLA") %>% 
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(index = ifelse(nobs >=20, 1,
                        ifelse(nobs < 20, 0.33+(nobs/30), 0)))

# this function for nobs < 20 has been calculated using the formula for a linear regression
# where index = a + b x nobs; 0,5 = a + bx5; 1 = a + bx20, which equals index = 0.33 + (1/30 x nobs)


## filter for Seed Mass, same process as above 
# replace those 0 values so they can be used
seed.now <- trait.cat %>% 
  filter(TraitShort == "SeedMass") %>% 
  mutate(StdValue = replace(StdValue, StdValue == 0, 0.0001)) %>%
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>% 
  distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(index = ifelse(nobs >=20, 1,
                      ifelse(nobs < 20, 0.33+(nobs/30), 0)))


## Additional step in here as we need to bring in the gap-filled data
gap.seeds <- read.csv("data/2022_gapfilled_seed_data.csv")

# extract unique columns so they can be combined
gap.unique1 <- gap.seeds %>% dplyr::select(sp, GapMedian, GapSD)

# extract the species vector with the gap-filled species
gap.vector1 <- unique(gap.seeds$sp)

# extract all the columns from the principal database above for the gap-filled species
gap.long1 <- trait.cat %>% 
  filter(sp %in% gap.vector1) %>% 
  dplyr::distinct(sp, .keep_all = TRUE)

# merge into one dataframe
gap.more1 <- merge(gap.long1, gap.unique1, by = "sp")
gap.more2 <- distinct(gap.more1, sp, .keep_all = TRUE)

# gap.seeds.new needs to have the same columns as seed.now so they can be binded
# we add an index value of 0.5 as we trust less the species that we have been gap-filling
# remove redundant columns to facilitate row bindng
# I am replacing the 0 values of the SD by half of the minimum value that is not zero
gap.seeds.new <- gap.more2 %>% 
  dplyr::mutate(LogTraitValue = NA) %>% 
  dplyr::mutate(TraitValueSD = NA) %>%
  dplyr::mutate(LogTraitValueSD = GapSD) %>% 
  dplyr::mutate(MedianTraitValue = NA) %>%
  dplyr::mutate(LogTraitMedianValue = GapMedian) %>% 
  dplyr::mutate(index = 0.5) %>%
  dplyr::mutate(TraitShort = "SeedMass") %>%
  dplyr::mutate(Dataset = "Gap-filled") %>% 
  dplyr::mutate(ValueKindName = "Gap-filled") %>% dplyr::mutate(DataContributor = "Gap-filled") %>% 
  dplyr::mutate(StdValue = NA) %>% dplyr::mutate(Source = "Gap-filled") %>%
  dplyr::select(., -c(GapMedian, GapSD))

# We have log-transformed seed values, calculated their SD and median, and then calculated the median of those
# for the gap-filling species. This is the closest possible way to do this similar to the species with data. 

# reorder columns in gap.seeds.new according to column order in seed.now
gap.seeds.new <- gap.seeds.new[names(seed.now)]

# bind dataframes into one: 40 obs (28 + 12)
seed.all <- bind_rows(seed.now, gap.seeds.new)
str(seed.all) 


  
## filter for Plant Height, same process as above
height.now <- trait.cat %>% 
  filter(TraitShort == "PlantHeight") %>% 
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>% 
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>% 
  distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(index = ifelse(nobs >=20, 1,
                      ifelse(nobs < 20, 0.33+(nobs/30), 0)))

hei.means <- height.now %>% group_by(gf) %>% summarise(mean(MedianTraitValue))



## SLA VARIATION MODEL ----
sla.var.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = sla.now, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_sla_var_mod")
summary(sla.var.mod) #negative ns
plot(sla.var.mod)


# Model predictions
newsladata2 = data.frame(LogTraitValueSD = seq(from = min(sla.now$LogTraitValueSD),
                                       to = max(sla.now$LogTraitValueSD), by = 0.1))
fit2 = fitted(
  sla.var.mod,
  newdata = newsladata2,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit2) = c('fit', 'se', 'lwr', 'upr')
df_plot2 = cbind(newsladata2, fit2)


# Black nice labels
(sla.var.plot <- ggplot() + 
    geom_point(data = sla.now, aes(x = LogTraitValueSD, y = CentredRangeLog, 
                                   colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("\nSLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_plot2, aes(x = LogTraitValueSD, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_plot2, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), 
                fill = "#0a6933", alpha = 0.5) +
    geom_label_repel(data = subset(sla.now, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                      "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitValueSD, CentredRangeLog, label = SpLabel), color = "black", box.padding = 2, 
                     segment.color = "black", fill = "white", label.size = 1,
                     fontface = "italic", size=8) +
    labs(colour = "Deciduousness") + range.theme + 
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))



## SEED VARIATION MODEL ----

# modeling including gap-filled species but fitting separate dispersal modes
berry.spp <- seed.all %>% filter(DispersalMode == "Berry")
wind.spp <- seed.all %>% filter(DispersalMode == "Wind")

# berry species
berry.seed.var.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = berry.spp, 
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_seed_var_berry_mod")
summary(berry.seed.var.mod) # negative ns

# wind species
wind.seed.var.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = wind.spp, 
                     iter = 2000, chains = 4, warmup = 400,
                     file = "models/2022_seed_var_wind_mod")
summary(wind.seed.var.mod) # positive ns


# Run the seed mass variation model without the gap-filled species (only true seed mass values)
seed.var.nogap.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = seed.now, 
                    iter = 2000, chains = 4, warmup = 400,
                    file = "models/2022_seed_var_nogap_mod")
summary(seed.var.nogap.mod) # negative ns



# Full model
seed.var.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = seed.all, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_seed_var_mod")

summary(seed.var.mod) # negative ns
#plot(seed.var.mod) 


## Model predictions
newseeddata2 = data.frame(LogTraitValueSD = seq(from = min(seed.all$LogTraitValueSD),
                                        to = max(seed.all$LogTraitValueSD), by = 0.1))
fitseed2 = fitted(
  seed.var.mod,
  newdata = newseeddata2,
  re_formula = NULL, # ignore random effects
  summary = TRUE) # mean and 95% CI

colnames(fitseed2) = c('fit', 'se', 'lwr', 'upr')
df_seedplot2 = cbind(newseeddata2, fitseed2)


# Let's filter depending on source
proper.sources <- c("LedaKleyer", "TTT", "TRY", "BobHollisterSeeds", "LucieSmrzova", "EstherLevesqueBylot", "No")
other.points <- filter(seed.all, Source %in% proper.sources)
gap.points <- filter(seed.all, Source == "Gap-filled")


# Plot seed mass & ranges
(plot.seed.now2 <- ggplot() + 
    geom_point(data = other.points,
               aes(x = LogTraitValueSD, y = CentredRangeLog, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points, shape = 1, 
               aes(x = LogTraitValueSD, y = CentredRangeLog, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    xlab("\nSeed Mass Variation (log mg)") + 
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seedplot2, aes(x = LogTraitValueSD, y = fit), colour = "#E57E00") +
    geom_ribbon(data = df_seedplot2, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), 
                fill = "#E57E00", alpha = 0.4) +
    geom_label_repel(data = subset(seed.all, sp %in% c("Dryas integrifolia", "Linnaea borealis", "Myrica gale")), 
                     aes(LogTraitValueSD, CentredRangeLog, label = SpLabel), 
                     color = "black", box.padding = 2, fill = "white",
                     segment.color = "black", label.size = 1,
                     fontface = "italic", size=8) +
    labs(colour = "Dispersal mode") + range.theme +
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))



## HEIGHT VARIATION MODEL ----
height.var.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitValueSD, data = height.now, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_height_var_mod")

summary(height.var.mod) #positive ns
plot(height.var.mod) #good convergence


## Model predictions
newheidata2 = data.frame(LogTraitValueSD = seq(from = min(height.now$LogTraitValueSD),
                                       to = max(height.now$LogTraitValueSD), by = 0.1))
fithei2 = fitted(
  height.var.mod,
  newdata = newheidata2,
  re_formula = NA, # ignore random effects
  summary = TRUE) # mean and 95% CI

colnames(fithei2) = c('fit', 'se', 'lwr', 'upr')
df_heiplot2 = cbind(newheidata2, fithei2)


## Plot height & ranges
(plot.height.now2 <- ggplot() +
    geom_point(data = height.now, aes(x = LogTraitValueSD, y = CentredRangeLog, 
                                      colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nHeight Variation (log m)") + 
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_line(data = df_heiplot2, aes(x = LogTraitValueSD, y = fit), colour = "#800080") +
    geom_ribbon(data = df_heiplot2, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), fill = "#800080", alpha = 0.3) +
    geom_label_repel(data = subset(height.now, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                         "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitValueSD, CentredRangeLog, label = SpLabel), fill = "white",
                     color = "black", box.padding = 2, segment.color = "black",
                     fontface = "italic", size=8, label.size = 1) +
    range.theme +
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
            legend.text=element_text(size = 20),
            legend.position = "top", legend.key = element_blank(),
            legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE)))





## FULL VARIATION MODEL WITH ALL TRAITS ---- 

# combine trait dataframes
all.traits <- rbind(sla.now, height.now, seed.all) 

# convert to long format
all.traits.long <- pivot_wider(all.traits, names_from = TraitShort, 
                               values_from = c(LogTraitValueSD, index))

# define function 
collapse <- function(x) x[!is.na(x)][1]

# merge rows per species and remove species without all 3 trait values (n = 34)
# calculate combined index by normalizing the individual trait indexes
all.traits.short <- all.traits.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()


# select the other relevant traits
more.traits <- all.traits %>% 
  dplyr::select(sp, CentredRangeLog, gf, Family, Genus, 
                LeafPhenology, DispersalMode, range_cat, category.abs) %>% 
  distinct(sp, .keep_all = TRUE)

# merge dataframes
full.traits <- left_join(all.traits.short, more.traits, by = "sp")
full.traits.yes <- distinct(full.traits, sp, .keep_all = TRUE)

# save this for the PCA (script #6)
write.csv(full.traits.yes, file = "data/2022_all_sp_trait_var.csv")



## Fit full model with all traits
all.traits.var.mod <- brm(CentredRangeLog | weights(combined.index) ~ 
                        Hei_var + Seed_var + SLA_var, data = full.traits.yes, 
                      iter = 2000, chains = 4, warmup = 400, 
                      file =  "models/2022_all_traits_var_mod")
summary(all.traits.var.mod) # ns
#plot(all.traits.var.mod) #model has converged well


## Fit full model with interactions on 2x2 
all.traits.var.mod3 <- brm(CentredRangeLog | weights(combined.index) ~ (Hei_var * Seed_var) + 
                             (Hei_var * SLA_var) + (SLA_var * Seed_var), data = full.traits.yes,
                           iter = 2000, chains = 4, warmup = 400, 
                           file = "models/2022_all_traits_var_int_mod")
summary(all.traits.var.mod3) # no significant interactions here
#plot(all.traits.var.mod3) #model has converged well




## SLA VALUES MODEL ----
sla.val.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, data = sla.now, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_sla_val_mod")
summary(sla.val.mod) # negative ns
#plot(sla.val.mod)


# model predictions
slaweidata2 <- data.frame(LogTraitMedianValue = 
                           seq(from = min(sla.now$LogTraitMedianValue),
                               to = max(sla.now$LogTraitMedianValue), by = 0.1))

fit.sla.wei2 = fitted(
  sla.val.mod,
  newdata = slaweidata2,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.sla.wei2) = c('fit', 'se', 'lwr', 'upr')
df_sla_wei2 = cbind(slaweidata2, fit.sla.wei2)


# Plot predictions
(predsla.wei <- ggplot() + 
    geom_point(data = sla.now, aes(x = LogTraitMedianValue, y = CentredRangeLog,
                                   colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("Median SLA value (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_sla_wei2, aes(x = LogTraitMedianValue, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_sla_wei2, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), 
                fill = "#0a6933", alpha = 0.5) +
    geom_label_repel(data = subset(sla.now, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                      "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitMedianValue, CentredRangeLog, label = SpLabel), color = "black", box.padding = 2, 
                     segment.color = "black", fill = "white", label.size = 1,
                     fontface = "italic", size=8) + range.theme +
    labs(colour = "Deciduousness") +
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))




## SEED VALUES MODEL ----

# modeling without the gap-filled species
seed.val.nogap.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, data = seed.now, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_seed_val_nogap_mod")
summary(seed.val.nogap.mod) # positive ns


# berry species
berry.seed.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, data = berry.spp, 
                      iter = 2000, chains = 4, warmup = 400,
                      file = "models/2022_seed_val_berry_mod")
summary(berry.seed.mod) # positive ns

# wind species
wind.seed.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, data = wind.spp, 
                      iter = 2000, chains = 4, warmup = 400,
                     file = "models/2022_seed_val_wind_mod")
summary(wind.seed.mod) # negative ns

# wind species have greater slopes that berry species, but not significantly.



# model including gap-filled species
seed.val.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, data = seed.all, 
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_seed_val_mod")
summary(seed.val.mod) # positive ns
#plot(seed.val.mod)


# model predictions
seedweidata <- data.frame(LogTraitMedianValue = 
                            seq(from = min(seed.all$LogTraitMedianValue),
                                to = max(seed.all$LogTraitMedianValue), by = 0.1))

fit.seed.wei = fitted(
  seed.val.mod,
  newdata = seedweidata,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.seed.wei) = c('fit', 'se', 'lwr', 'upr')
df_seed_wei = cbind(seedweidata, fit.seed.wei)


# Plot predictions
(predseed.wei <- ggplot() + 
    geom_point(data = other.points, 
               aes(x = LogTraitMedianValue, y = CentredRangeLog, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points, shape = 1, 
               aes(x = LogTraitMedianValue, y = CentredRangeLog, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    xlab("\nMedian Seed Mass value (log mg)") + 
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seed_wei, aes(x = LogTraitMedianValue, y = fit), colour = "#E57E00") + 
    geom_ribbon(data = df_seed_wei, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), 
                fill = "#E57E00", alpha = 0.4) + 
    labs(colour = "Dispersal mode") +
    geom_label_repel(data = subset(seed.all, sp %in% c("Dryas integrifolia", "Linnaea borealis", "Myrica gale")), 
                     aes(LogTraitMedianValue, CentredRangeLog, label = SpLabel), 
                     color = "black", box.padding = 2, fill = "white",
                     segment.color = "black", label.size = 1,
                     fontface = "italic", size=8) + range.theme +
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))




## HEIGHT VALUES MODEL ----
hei.val.mod <- brm(CentredRangeLog | weights(index) ~ LogTraitMedianValue, 
                   data = height.now, iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_hei_val_mod")
summary(hei.val.mod) # negative ns


# model predictions
heiweidata <- data.frame(LogTraitMedianValue = 
                            seq(from = min(height.now$LogTraitMedianValue),
                                to = max(height.now$LogTraitMedianValue), by = 0.1))
fit.hei.wei = fitted(
  hei.val.mod,
  newdata = heiweidata,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.hei.wei) = c('fit', 'se', 'lwr', 'upr')
df_hei_wei = cbind(heiweidata, fit.hei.wei)


# Plot predictions
(predhei.wei <- ggplot() + 
    geom_point(data = height.now, aes(x = LogTraitMedianValue, y = CentredRangeLog, 
                                      colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nMedian Height value (log m)") + 
    ylab(expression(bold(paste("Current Species Range (log km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_line(data = df_hei_wei, aes(x = LogTraitMedianValue, y = fit), colour = "#800080") + 
    geom_ribbon(data = df_hei_wei, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), 
                fill = "#800080", alpha = 0.3) +
    geom_label_repel(data = subset(height.now, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                         "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitMedianValue, CentredRangeLog, label = SpLabel), fill = "white",
                     color = "black", box.padding = 2, segment.color = "black",
                     fontface = "italic", size=8, label.size = 1) + range.theme +
    theme(axis.ticks.length = unit(.25, "cm"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE)))




## FULL VALUES MODEL ---- 

# convert to long format both the value and the index
all.traits.val <- pivot_wider(all.traits, 
                               names_from = TraitShort, 
                               values_from = c(LogTraitMedianValue, index))

# combine the trait values and calculate combined index
all.traits.val.short <- all.traits.val %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   SeedMass_val = collapse(LogTraitMedianValue_SeedMass), 
                   SLA_val = collapse(LogTraitMedianValue_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# merge dataframes
full.traits.val <- left_join(all.traits.val.short, more.traits, by = "sp")
full.traits.val.yes <- distinct(full.traits.val)

# save this for the PCA (script #6)
write.csv(full.traits.val.yes, file = "data/2022_all_sp_trait_val.csv")


# fit model with weights() function and combined index
all.traits.val.mod <- brm(CentredRangeLog | weights(combined.index) ~ 
                            SLA_val + Hei_val + SeedMass_val, 
                          data = full.traits.val.yes, iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_all_traits_val_mod")

summary(all.traits.val.mod) # all ns
#plot(all.traits.val.mod)



# model with 2x2 interactions 
all.traits.val.mod4 <- brm(CentredRangeLog | weights(combined.index) ~ 
                             (SLA_val * Hei_val) + (SeedMass_val * SLA_val) + (SeedMass_val * Hei_val), 
                           data = full.traits.val.yes, 
                           iter = 2000, chains = 4, warmup = 400, 
                           file = "models/2022_all_traits_val_int_mod")
summary(all.traits.val.mod4)  # all ns



## FIGURE 3 (FULL PANEL) ----
(current.panel <- ggarrange(predhei.wei, predsla.wei, predseed.wei,
                            plot.height.now2, sla.var.plot, plot.seed.now2,
                            nrow = 2, ncol = 3, 
                            labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                            font.label = list(size = 30)))

ggplot2::ggsave(current.panel, filename = "figures/Figure_3.jpeg", 
                width = 60, height = 60, units = "cm", dpi = 500)




## TRAIT CORRELATIONS ----

hei.sla <- brm(SLA_val ~ Hei_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400, 
               file = "models/2022_sla_hei_cor_mod")
summary(hei.sla) # positive ns

sla.hei <- brm(Hei_val ~ SLA_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400,
               file = "models/2022_hei_sla_cor_mod")
summary(hei.sla) # positive ns


hei.seed <- brm(SeedMass_val ~ Hei_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400,
                file = "models/2022_hei_seed_cor_mod")
summary(hei.seed) # positive ns

seed.hei <- brm(Hei_val ~ SeedMass_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400,
                file = "models/2022_seed_hei_cor_mod")
summary(seed.hei) # positive ns


seed.sla <- brm(SeedMass_val ~ SLA_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400,
                file = "models/2022_seed_sla_cor_mod")
summary(seed.sla) # positive ns

sla.seed <- brm(SLA_val ~ SeedMass_val, data = all.traits.val.short, iter = 2000, chains = 4, warmup = 400,
                file = "models/2022_sla_seed_cor_mod")
summary(sla.seed) # negative ns





## CATEGORICAL MODELS ----


## SLA variation
sla.var.cat.mod <- brm(LogTraitValueSD | weights(index) ~ category.abs, data = sla.now, 
                       iter = 2000, chains = 4, warmup = 400,
                       file = "models/2022_sla_var_cat_mod")
summary(sla.var.cat.mod) # no significant difference between the three categories
conditional_effects(sla.var.cat.mod)

## Seed variation
seed.var.cat.mod <- brm(LogTraitValueSD | weights(index) ~ category.abs, data = seed.all, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_seed_var_cat_mod")
summary(seed.var.cat.mod) # no significant difference between the three categories
conditional_effects(seed.var.cat.mod)


## Height variation
hei.var.cat.mod <- brm(LogTraitValueSD | weights(index) ~ category.abs, data = height.now, 
                        iter = 2000, chains = 4, warmup = 400,
                       file = "models/2022_hei_var_cat_mod")
summary(hei.var.cat.mod) # no significant difference between the three categories
conditional_effects(hei.var.cat.mod)


## SLA values
sla.val.cat.mod <- brm(LogTraitMedianValue | weights(index) ~ category.abs, data = sla.now, 
                       iter = 2000, chains = 4, warmup = 400,
                       file = "models/2022_sla_val_cat_mod")
summary(sla.val.cat.mod) # no significant difference between the three categories
conditional_effects(sla.val.cat.mod)


## Seed values
seed.val.cat.mod <- brm(LogTraitMedianValue | weights(index) ~ category.abs, data = seed.all, 
                       iter = 2000, chains = 4, warmup = 400,
                       file = "models/2022_seed_val_cat_mod")
summary(seed.val.cat.mod) # no significant difference between the three categories
conditional_effects(seed.val.cat.mod)


## Height values
hei.val.cat.mod <- brm(LogTraitMedianValue | weights(index) ~ category.abs, data = height.now, 
                        iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_hei_val_cat_mod")
summary(hei.val.cat.mod) # no significant difference between the three categories
conditional_effects(hei.val.cat.mod)





# For the models below we need all the species, not only those that have all 3 trait values (traits don't matter here)
cur.cat <- all.traits %>% select(CentredRangeLog, sp, gf, DispersalMode, LeafPhenology, Family) %>% 
  distinct(sp, .keep_all = TRUE) #62 spps

## Current ranges & functional group
current.fg.mod <- brm(CentredRangeLog ~ gf, data = cur.cat, 
                      iter = 2000, chains = 4, warmup = 400,
                      file = "models/2022_cur_fg_mod")
summary(current.fg.mod) # no significant difference between the three functional groups
conditional_effects(current.fg.mod)


## Current ranges & dispersal mode
current.disp.mod <- brm(CentredRangeLog ~ DispersalMode, data = cur.cat, 
                        iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_cur_disp_mod")
summary(current.disp.mod) # no significant difference between the three dispersal modes
conditional_effects(current.disp.mod) # but it's interesting that wind-dispersed have smaller current ranges


## Current ranges & decidiousness
current.decid.mod <- brm(CentredRangeLog ~ LeafPhenology, data = cur.cat, 
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_cur_decid_mod")
summary(current.decid.mod) # no significant difference between evergreen and deciduous
conditional_effects(current.decid.mod)


## Current ranges & family
current.fam.mod <- brm(CentredRangeLog ~ Family, data = cur.cat, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_cur_fam_mod")
summary(current.fam.mod) # no significant difference between the 7 families
conditional_effects(current.fam.mod) # Except for Salicaceae and Rosaceae, which overlap with all other families
# except with each other. Salicaceae has smaller current ranges than Rosaceae.



## BINOMIAL MODELS ----

## Category as a function of trait variation - binomial and weights

# Loser and no changers go together in the same category (40 obs)
full.traits.var.yes.bi <- full.traits.yes %>% 
  mutate(category.abs.nb = ifelse(category.abs == "Loser", 0, 
                                  ifelse(category.abs == "No change", 0, 1)))


all.traits.var.mod25 <- brm(category.abs.nb | weights(combined.index) ~ SLA_var + Hei_var + Seed_var, 
                           data = full.traits.var.yes.bi, iter = 2000, chains = 4, 
                           warmup = 400, family = bernoulli(link = "logit"),
                           file = "models/2022_all_traits_var_cat")
summary(all.traits.var.mod25) # seed variation positive significant
plot(conditional_effects(all.traits.var.mod25), points = TRUE)



# extract fixed effects
coefs.var.fit <- as.data.frame(fixef(all.traits.var.mod25))
coefs.var.fit2 <- coefs.var.fit %>% rownames_to_column("Traits")
coefs.var.fit3 <- coefs.var.fit2[-1,]

coefs.var.fit3$Traits[coefs.var.fit3$Traits == 'SLA_var'] <- 'SLA'
coefs.var.fit3$Traits[coefs.var.fit3$Traits == 'Hei_var'] <- 'Plant Height'
coefs.var.fit3$Traits[coefs.var.fit3$Traits == 'Seed_var'] <- 'Seed Mass'


# Plot posteriors and CIs
(coefs.var.plot <- ggplot(coefs.var.fit3, aes(x = Traits, y = Estimate, colour = Traits)) + 
    geom_point(size = 10) + 
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), size = 1.5, width = 0.5) + 
    ylab("Binomial model effect sizes\nCategory vs trait variation\n") + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    annotate(geom="text", x = 2, y = 4, label="*", size = 15) +
    scale_colour_manual(values = c("#800080", "#E57E00", "#0a6933")) + range.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text=element_text(size = 16), plot.margin = unit(c(2,2,2,2), "cm")))



## Category as a function of trait values - binomial and weights

# Loser and no changers go together in the same category
full.traits.val.yes.bi.23 <- full.traits.val.yes %>% 
  mutate(category.abs.nb = ifelse(category.abs == "Loser", 0, 
                                  ifelse(category.abs == "No change", 0, 1)))


all.traits.val.mod234 <- brm(category.abs.nb | weights(combined.index) ~ SLA_val + Hei_val + SeedMass_val, 
                            data = full.traits.val.yes.bi.23, iter = 2000, chains = 4, 
                            warmup = 400, family = bernoulli(link = "logit"), 
                            file = "models/2022_all_traits_val_cat")
summary(all.traits.val.mod234) #all ns
plot(all.traits.val.mod234)
plot(conditional_effects(all.traits.val.mod234), points = TRUE)



## Extract fixed effects
coefs.val.fit <- as.data.frame(fixef(all.traits.val.mod234))
coefs.val.fit2 <- coefs.val.fit %>% rownames_to_column("Traits")
coefs.val.fit3 <- coefs.val.fit2[-1,]

coefs.val.fit3$Traits[coefs.val.fit3$Traits == 'SLA_val'] <- 'SLA'
coefs.val.fit3$Traits[coefs.val.fit3$Traits == 'Hei_val'] <- 'Plant Height'
coefs.val.fit3$Traits[coefs.val.fit3$Traits == 'SeedMass_val'] <- 'Seed Mass'


# Plot posteriors and CIs
(coefs.val.plot <- ggplot(coefs.val.fit3, aes(x = Traits, y = Estimate, colour = Traits)) + 
    geom_point(size = 10) + 
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), size = 1.5, width = 0.5) + 
    ylab("Binomial model effect sizes\nCategory vs trait values\n") + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_colour_manual(values = c("#800080", "#E57E00", "#0a6933")) + range.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text=element_text(size = 16), plot.margin = unit(c(2,2,2,2), "cm")))
    







## BINOMIAL MODELS FOR COVER CHANGE CATEGORY ----

# Loser and no changers go together in the same category
cov.change.cat <- read.csv("data/CoverChangeSpecies_QuantsCats_Final.csv")
names(cov.change.cat)[names(cov.change.cat) == 'Name'] <- 'sp' # (36 spps with cover data)


## Trait variation (30 spp with cover data and traits)
full.traits.var.cov.bi <- left_join(full.traits.yes, cov.change.cat, by = "sp")
full.traits.var.cov.bi2 <- full.traits.var.cov.bi %>% drop_na(CoverCategory) %>%
  mutate(category.cov = ifelse(CoverCategory == "Loser", 0, 
                                  ifelse(CoverCategory == "No change", 0, 1)))


# binomial model with loser and no change in the same group and winner in the other
all.traits.var.cov.mod <- brm(category.cov | weights(combined.index) ~ SLA_var + Hei_var + Seed_var, 
                            data = full.traits.var.cov.bi2, iter = 2000, chains = 4, 
                            warmup = 400, family = bernoulli(link = "logit"),
                            file = "models/2022_all_traits_var_cat_cov")
summary(all.traits.var.cov.mod) # all ns
plot(conditional_effects(all.traits.var.cov.mod), points = TRUE)



## Trait values
full.traits.val.cov.bi <- left_join(full.traits.val.yes, cov.change.cat, by = "sp")
full.traits.val.cov.bi2 <- full.traits.val.cov.bi %>% drop_na(CoverCategory) %>%
  mutate(category.cov = ifelse(CoverCategory == "Loser", 0, 
                               ifelse(CoverCategory == "No change", 0, 1)))


# Binomial model 
all.traits.val.cov.mod <- brm(category.cov | weights(combined.index) ~ SLA_val + Hei_val + SeedMass_val, 
                             data = full.traits.val.cov.bi2, iter = 2000, chains = 4, 
                             warmup = 400, family = bernoulli(link = "logit"), 
                             file = "models/2022_all_traits_val_cat_cov")
summary(all.traits.val.cov.mod) # all ns
plot(all.traits.val.cov.mod)
plot(conditional_effects(all.traits.val.cov.mod), points = TRUE)



## COVER CHANGE FIGURE ----

# Slope variation per subsite
cov.change.slopes <- read.csv("data/CoverChange_AllSpecies_NewModel_Slopes.csv")
str(cov.change.slopes)
colnames(cov.change.slopes)[which(names(cov.change.slopes) == "X5.")] <- "Quantile5"
colnames(cov.change.slopes)[which(names(cov.change.slopes) == "X50.")] <- "Quantile50"
colnames(cov.change.slopes)[which(names(cov.change.slopes) == "X95.")] <- "Quantile95"
colnames(cov.change.slopes)[which(names(cov.change.slopes) == "Name")] <- "sp"

cov.change.slopes$sp <- as.character(cov.change.slopes$sp)
cov.change.slopes$sp[cov.change.slopes$sp == "Ledum palustre"] <- "Rhododendron tomentosum"
cov.change.slopes <- filter(cov.change.slopes, X != 132)
cov.change.cat$sp <- as.character(cov.change.cat$sp)
cov.change.cat$sp[cov.change.cat$sp == "Ledum palustre"] <- "Rhododendron tomentosum"


# Order in descending cover change order
win.to.lose <- cov.change.cat %>% arrange(desc(CoverCategory), desc(MeanSlope))
win.to.lose.vector <- unique(win.to.lose$sp)


# Here comes the plot (Fig S3)
(cov.slopes.plot <-  ggplot() +
    geom_point(data = cov.change.slopes, aes(x = factor(sp, level = win.to.lose.vector), y = Quantile50), 
                                             colour = "lightgrey", size = 3, alpha = 0.7) + 
    geom_errorbar(data = cov.change.slopes, aes(x = factor(sp, level = win.to.lose.vector), ymin = Quantile5, ymax = Quantile95), 
                  colour = "lightgrey", width = 0.2, size = 0.5) +
    geom_point(data = cov.change.cat, aes(x = factor(sp, level = win.to.lose.vector), y = MeanSlope, colour = CoverCategory), 
                                          size = 5.5) + 
    geom_errorbar(data = cov.change.cat, aes(x = factor(sp, level = win.to.lose.vector), ymin = LowerQuantile, ymax = UpperQuantile), 
                                             colour = "black", size = 0.5, width = 0.7) + 
    scale_colour_manual(values = c("#CC6677", "#9eabad", "#44AA99")) + geom_hline(yintercept = 0) +
    xlab("\nSpecies") + ylab("Cover change over time\nEffect sizes\n") + 
    range.theme + labs(colour = "Category") +
    theme(axis.title.x = element_text(face="bold", size=22),
          axis.text.x  = element_text(face = "italic", size=15, colour = "black", angle = 52, vjust = 1, hjust = 1),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))

ggsave(cov.slopes.plot, filename = "figures/Figure_S3.png", 
       width = 30, height = 20, units = "cm")
