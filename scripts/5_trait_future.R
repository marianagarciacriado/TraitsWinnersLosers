## Trait-range manuscript
## Mariana Garcia
## March 2019
## Script 5. Relationships between traits and species range changes


## LIBRARIES ----
library(tidyverse)
library(ggpubr)
library(rstan)
library(brms)
library(ggrepel)
library(dplyr)


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



## DATA PREP ----
load("data/2022_trait_cat.RData")
#this is assigned to the object 'trait.cat' - 17921 obs


## Create objects by filtering per trait

# filter for SLA and log-transform the individual trait records
# group per species, calculate SD and median 
sla.fut <- trait.cat %>% 
  dplyr::filter(TraitShort == "SLA") %>% 
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>%
  dplyr::distinct(sp, .keep_all = TRUE) %>% 
  ungroup() %>%
  mutate(index = ifelse(nobs >=20, 1,
                        ifelse(nobs < 20, 0.33+(nobs/30), 0)))

# this function for nobs < 30 has been calculated using the formulae for a linear regession
# where index = a + b x nobs; 0,5 = a + bx5; 1 = a + bx20, which equals index = 0.33 + (1/30 x nobs)

# check all distributions of values for starters
hist(sla.fut$rel.median) # normal distrib kind of around zero but with long tail
hist(sla.fut$rel.median.log) # closer to zero and kind of around it
hist(sla.fut$cent.rel.median.log) # values are now centered around zero

hist(sla.fut$abs.median) # normal distrib kind of around zero but with long tail
hist(sla.fut$abs.median.log) # two tails at both sides of zero
hist(sla.fut$cent.abs.median.log) # two tails at both sides of zero

# check SLA standard deviation
hist(sla.fut$TraitValueSD) #normal distribution, not too far from zero
hist(sla.fut$LogTraitValueSD) #much closer to zero

# check SLA values
hist(sla.fut$MedianTraitValue) #normal distribution, not too far from zero
hist(sla.fut$LogTraitMedianValue) #much closer to zero



## Plant height
hei.fut <- trait.cat %>% 
  dplyr::filter(TraitShort == "PlantHeight") %>% 
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>% 
  dplyr::distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>% 
  mutate(index = ifelse(nobs >=20, 1,
                      ifelse(nobs < 20, 0.33+(nobs/30), 0)))

# check all distributions of values for starters
hist(hei.fut$rel.median) # normal distrib kind of around zero but with long tail
hist(hei.fut$rel.median.log) # closer to zero but still not around zero
hist(hei.fut$cent.rel.median.log) # values are now centered around zero

hist(hei.fut$abs.median) # normal distrib kind of around zero but with long tail
hist(hei.fut$abs.median.log) # closer to zero but still not around zero
hist(hei.fut$cent.abs.median.log) # values are now centered around zero

# check Height standard deviation
hist(hei.fut$TraitValueSD) #normal distribution, not too far from zero
hist(hei.fut$LogTraitValueSD) #much closer to zero

# check Height values
hist(hei.fut$MedianTraitValue) #normal distribution, not too far from zero
hist(hei.fut$LogTraitMedianValue) #much closer to zero



## Seed Mass
seed.fut <- trait.cat %>% 
  dplyr::filter(TraitShort == "SeedMass") %>% 
  mutate(StdValue = replace(StdValue, StdValue == 0, 0.0001)) %>%
  mutate(LogTraitValue = log(StdValue)) %>%
  group_by(sp) %>% 
  mutate(TraitValueSD = sd(StdValue)) %>%
  mutate(LogTraitValueSD = sd(LogTraitValue)) %>%
  mutate(MedianTraitValue = median(StdValue)) %>%
  mutate(LogTraitMedianValue = median(LogTraitValue)) %>% 
  dplyr::distinct(sp, .keep_all = TRUE) %>%
  ungroup() %>% 
  mutate(index = ifelse(nobs >=20, 1,
                        ifelse(nobs < 20, 0.33 +(nobs/30), 0)))


## Additional step in here as we need to bring in the gap-filled data
gap.seeds.fut <- read.csv("data/2022_gapfilled_seed_data.csv")

# extract unique columns so they can be combined
gap.unique <- gap.seeds.fut %>% dplyr::select(sp, GapMedian, GapSD)

# extract the species vector with the gap-filled species
gap.vector <- unique(gap.seeds.fut$sp)

# extract all the range columns from the principal database above for the gap-filled species
gap.long <- trait.cat %>% filter(sp %in% gap.vector) %>% distinct(sp, .keep_all = TRUE)

# merge into one dataframe
gap.more <- merge(gap.long, gap.unique, by = "sp")

# should be 12 obs
gap.more.more <- distinct(gap.more, sp, .keep_all = TRUE)

# gap.seeds.fut needs to have the same columns as seed.fut so they can be binded
# we add an index value of 0.5 as we trust less the species that we have been gap-filling
# remove redundant columns to facilitate row bindng
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

# We have log-transformed seed values, calculated their SD and median, and then calculated the median of those
# for the gap-filling species. This is the closest possible way to do this similar to the species with data. 

# reorder columns in gap.seeds.new according to column order in seed.now
gap.seeds.fut.new2 <- gap.seeds.fut.new[names(seed.fut)]

# bind dataframes into one
seed.fut.all <- bind_rows(seed.fut, gap.seeds.fut.new2)
str(seed.fut.all)




## SLA VARIATION / ABSOLUTE CHANGE ----

## SLA - 25% quantile model
sla.abs25.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitValueSD, data = sla.fut, 
                     iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_25abs_sla_var_mod")
summary(sla.abs25.mod) # positive ns
plot(sla.abs25.mod)



## SLA - median model
sla.medabs.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = sla.fut,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_medabs_sla_var_mod")
summary(sla.medabs.mod) # positive ns
plot(sla.medabs.mod) 


## Model predictions
sla.abs.data = data.frame(LogTraitValueSD = seq(from = min(sla.fut$LogTraitValueSD),
                                                to = max(sla.fut$LogTraitValueSD), by = 0.1))

fit.sla.med.abs = fitted(
  sla.medabs.mod,
  newdata = sla.abs.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.sla.med.abs) = c('fit', 'se', 'lwr', 'upr')
df_slamedabs = cbind(sla.abs.data, fit.sla.med.abs)


## Plot relationships
(sla.absmed.plot <- ggplot() + 
    geom_point(data = sla.fut, aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                                    colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("\nSLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_slamedabs, aes(x = LogTraitValueSD, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_slamedabs, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), 
                fill = "#0a6933", alpha = 0.5) +
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    geom_label_repel(data = subset(sla.fut, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                      "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitValueSD, cent.abs.median.log, label = SpLabel), color = "black", box.padding = 2, 
                     segment.color = "black", fill = "white", label.size = 1,
                     fontface = "italic", size=8) +
    range.theme + labs(colour = "Deciduousness") +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))


## SLA - 75% quantile model
sla.abs75.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitValueSD, 
                     data = sla.fut, iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_75abs_sla_var_mod")
summary(sla.abs75.mod)
plot(sla.abs75.mod) # positive significant


## Model predictions
sla.abs.data2 = data.frame(LogTraitValueSD = seq(from = min(sla.fut$LogTraitValueSD),
                                                to = max(sla.fut$LogTraitValueSD), by = 0.05))

fit.sla.abs75 = fitted(
  sla.abs75.mod,
  newdata = sla.abs.data2,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.sla.abs75) = c('fit', 'se', 'lwr', 'upr')
df_slaabs75 = cbind(sla.abs.data2, fit.sla.abs75)


## Back-transforming for plotting

# Mean value used for centering data back in file range.full2
range.full2 <- read.csv("data/range_full2.csv")
mean(range.full2$abs.quan75.log) #2.632218

# Mean value of the predicted values for Y (fit)
mean(df_slaabs75$fit) # 0.09216923

# Back-centre and back-log-transform
df_slaabs75_more <- df_slaabs75 %>% mutate(y_constant = 2.632218+fit) %>% 
  mutate(y_constant_exp = exp(y_constant)) %>% mutate(back_x = exp(LogTraitValueSD))

# Mean value of the un-centred Y values
mean(df_slaabs75_more$y_constant) # 2.724387

# Mean value of the un-centred un-logged Y values
mean(df_slaabs75_more$y_constant_exp) # 15.43348

# Back-transformed mean value of the un-centred un-logged Y values
exp(2.724387) # 15.24706

# Min abs + 1 for absolute values = 11.531 +1 = 12.531
df_slaabs75_more_more <- df_slaabs75_more %>% 
  mutate(y_constant_exp_offset = y_constant_exp - 12.531) %>% 
  mutate(y_cnt_exp_offset_mil = y_constant_exp_offset*1000000)

# Fix the issue with the back-transformation here
sla.fut2 <- sla.fut %>% mutate(ExpTraitValueSD = exp(LogTraitValueSD))

# Separate points per weighting
down.sla1 <- filter(sla.fut2, nobs <20)
up.sla1 <- filter(sla.fut2, nobs >=20)

# Plot this
(back.sla.abs75.var.plot2 <- ggplot() + 
    geom_point(data = df_slaabs75_more_more, aes(x = back_x, y = y_cnt_exp_offset_mil), colour = "black", size = 4) +
    geom_point(data = down.sla1, aes(x = ExpTraitValueSD, y = AbsoluteRangeChangeKm), colour = "#0a6933", size = 4, alpha = 0.7) + 
    geom_point(data = up.sla1, aes(x = ExpTraitValueSD, y = AbsoluteRangeChangeKm), colour = "#0a6933", size = 7, alpha = 0.7) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    xlab(expression(bold(paste("SLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("75% Absolute Species Range Change (km\n"^bold("2\n"), ")\n")))) + range.theme +
    theme(axis.title.x = element_text(face="bold", size=22),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.ticks.length = unit(.25, "cm")))




## HEIGHT VARIATION / ABSOLUTE CHANGE ----

## Height - 25% quantile model
hei.abs25.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitValueSD, data = hei.fut,
                     iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_25abs_hei_var_mod")
summary(hei.abs25.mod) # positive ns
plot(hei.abs25.mod)


## Height - median model 
hei.medabs.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = hei.fut,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_medabs_hei_var_mod")
summary(hei.medabs.mod) # positive ns
plot(hei.medabs.mod)

# Model predictions
hei.abs.data = data.frame(LogTraitValueSD = seq(from = min(hei.fut$LogTraitValueSD),
                                                to = max(hei.fut$LogTraitValueSD), by = 0.1))

fit.heimedabs = fitted(
  hei.medabs.mod,
  newdata = hei.abs.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.heimedabs) = c('fit', 'se', 'lwr', 'upr')
df_heimedabs = cbind(hei.abs.data, fit.heimedabs)

# Plot relationships
(hei.medabs.plot <- ggplot() + 
    geom_point(data = hei.fut, aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                                   colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nHeight variation (log m)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    geom_line(data = df_heimedabs, aes(x = LogTraitValueSD, y = fit), colour = "#800080") + 
    geom_ribbon(data = df_heimedabs, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr), fill = "#800080", alpha = 0.3) +
    range.theme +
    geom_label_repel(data = subset(hei.fut, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                         "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitValueSD, cent.abs.median.log, label = SpLabel), fill = "white",
                     color = "black", box.padding = 2, segment.color = "black",
                     fontface = "italic", size=8, label.size = 1) +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2, nrow=2, byrow=TRUE)))



## Height - 75% quantile model
hei.abs75.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitValueSD, 
                     data = hei.fut, iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_75abs_hei_var_mod")
summary(hei.abs75.mod) # positive ns
plot(hei.abs75.mod) 




## SEED MASS VARIATION / ABSOLUTE CHANGE ----

## Seed mass - 25% quantile model
seed.abs25.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitValueSD, 
                      data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_25abs_seed_var_mod")
print(seed.abs25.mod, digits = 5) # positive significant
plot(seed.abs25.mod)


## Dispersal mode with range shifts

# modeling including gap-filled species but fitting separate dispersal modes
berry.spp.rc <- seed.fut.all %>% filter(DispersalMode == "Berry")
wind.spp.rc <- seed.fut.all %>% filter(DispersalMode == "Wind")

# berry species
berry.seed.var.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = berry.spp.rc, 
                          iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_seed_var_berry_rc_mod")
summary(berry.seed.var.rc.mod) # negative ns

# wind species
wind.seed.var.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = wind.spp.rc, 
                         iter = 2000, chains = 4, warmup = 400,
                         file = "models/2022_seed_var_wind_rc_mod")
summary(wind.seed.var.rc.mod) # positive significant 


# Run the seed mass variation model without the gap-filled species (only true seed mass values)
seed.var.nogap.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, data = seed.fut, 
                          iter = 2000, chains = 4, warmup = 400,
                          file = "models/2022_seed_var_nogap_rc_mod")
summary(seed.var.nogap.rc.mod) # positive ns




## Seed median model
seed.medabs.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD,
                       data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_medabs_seed_var_mod")
summary(seed.medabs.mod) # positive significant
plot(seed.medabs.mod) 


# Model predictions
seed.abs.data = data.frame(LogTraitValueSD = seq(from = min(seed.fut.all$LogTraitValueSD),
                                                 to = max(seed.fut.all$LogTraitValueSD), by = 0.1))

fit.seedmedabs = fitted(
  seed.medabs.mod,
  newdata = seed.abs.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.seedmedabs) = c('fit', 'se', 'lwr', 'upr')
df_seedmedabs = cbind(seed.abs.data, fit.seedmedabs)


# Let's filter depending on source
proper.sources <- c("LedaKleyer", "TTT", "TRY", "BobHollisterSeeds", "LucieSmrzova", "EstherLevesqueBylot", "No")
other.points.fut <- filter(seed.fut.all, Source %in% proper.sources)
gap.points.fut <- filter(seed.fut.all, Source == "Gap-filled")



# Plot relationships
(seed.medabs.plot <- ggplot() + 
    geom_point(data = other.points.fut,
               aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points.fut, shape = 1, 
               aes(x = LogTraitValueSD, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    xlab("\nSeed Mass Variation (log mg)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seedmedabs, aes(x = LogTraitValueSD, y = fit), colour = "#E57E00") + 
    geom_ribbon(data = df_seedmedabs, aes(x = LogTraitValueSD, ymin = lwr, ymax = upr),
                fill = "#E57E00", alpha = 0.4) +
    geom_label_repel(data = subset(seed.fut.all, sp %in% c("Dryas integrifolia", "Myrica gale", "Linnaea borealis")), 
                     aes(LogTraitValueSD, cent.abs.median.log, label = SpLabel), 
                     color = "black", box.padding = 2, fill = "white",
                     segment.color = "black", label.size = 1,
                     fontface = "italic", size=8) +
    labs(colour = "Dispersal mode") + range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))



## Seed 75% quantile model
seed.abs75.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitValueSD,
                      data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_75abs_seed_var_mod")

summary(seed.abs75.mod) # positive significant



## FULL MODEL - ABSOLUTE/VARIATION ----

# combine trait dataframes
all.traits4 <- rbind(sla.fut, hei.fut, seed.fut.all)

# convert to long format
all.traits.long4 <- tidyr::pivot_wider(all.traits4, names_from = TraitShort, 
                                values_from = c(LogTraitValueSD, index))

# define function 
collapse <- function(x) x[!is.na(x)][1]

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.short4 <- all.traits.long4 %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# select the other relevant traits
more.traits4 <- all.traits4 %>% 
  dplyr::select(sp, cent.abs.median.log, cent.rel.median.log, gf, 
                Family, Genus, LeafPhenology, DispersalMode, category.abs, range_cat)

# merge dataframes
full.traits4 <- left_join(all.traits.short4, more.traits4, by = "sp")
full.traits.yes4 <- distinct(full.traits4, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.var.abs.mod <- brm(cent.abs.median.log | weights(combined.index) ~ Hei_var + Seed_var + SLA_var,
                              data = full.traits.yes4, iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_abs_var_mod")

summary(all.traits.var.abs.mod) # no variable is significant
plot(all.traits.var.abs.mod)


## Fit full model with all traits and 2x2 interaction 
all.traits.var.abs.mod4 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                 (Hei_var * Seed_var) + (SLA_var * Hei_var) + (SLA_var * Seed_var),
                               data = full.traits.yes4, iter = 2000, chains = 4, warmup = 400, 
                               file = "models/2022_all_traits_int_abs_var_mod")

summary(all.traits.var.abs.mod4) # no variable is significant
conditional_effects(all.traits.var.abs.mod4)



## SLA VARIATION / RELATIVE CHANGE ----

## SLA - 25% quantile model
sla.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitValueSD, data = sla.fut, 
                      iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_25rel_sla_var_mod")
summary(sla.rel25.mod) # positive significant
plot(sla.rel25.mod) 



## SLA - median model
sla.medrel.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, data = sla.fut, 
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_medrel_sla_var_mod")
summary(sla.medrel.mod) # positive significant
plot(sla.medrel.mod) 



## SLA - 75% quantile model
sla.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitValueSD, data = sla.fut, 
                         iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_75rel_sla_var_mod")
summary(sla.rel75.mod) # positive significant 
plot(sla.rel75.mod) 


# Model predictions
sla.rel75.data = data.frame(LogTraitValueSD = seq(from = min(sla.fut$LogTraitValueSD),
                                                 to = max(sla.fut$LogTraitValueSD), by = 0.05))

fit.slarel75 = fitted(
  sla.rel75.mod,
  newdata = sla.rel75.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.slarel75) = c('fit', 'se', 'lwr', 'upr')
df_slarel75 = cbind(sla.rel75.data, fit.slarel75)


## Back-transforming 

# Mean value used for centering data back in file range.full2
mean(range.full2$rel.quan75.log) #0.7739364

# Mean value of the predicted values for Y (fit)
mean(df_slarel75$fit) # 0.114976

# Back-centre and back-log-transform
df_slarel75_more <- df_slarel75 %>% mutate(y_constant = 0.7739364+fit) %>% 
  mutate(y_constant_exp = exp(y_constant)) %>% mutate(back_x = exp(LogTraitValueSD))

# Mean value of the un-centred Y values
mean(df_slarel75_more$y_constant) # 0.8889124

# Mean value of the un-centred un-logged Y values
mean(df_slarel75_more$y_constant_exp) # 2.472609

# Back-transformed mean value of the un-centred un-logged Y values
exp(2.472609) # 11.85333

# Min abs + 1 for absolute values = 1+1 = 2
df_slarel75_more_more <- df_slarel75_more %>% 
  mutate(y_constant_exp_offset = y_constant_exp - 2) %>% 
  mutate(y_cnt_exp_offset_mil = y_constant_exp_offset*100)

# Plot 
(back.sla.rel75.var.plot2 <- ggplot() + 
    geom_point(data = df_slarel75_more_more, aes(x = back_x, y = y_cnt_exp_offset_mil), colour = "black", size = 4) +
    geom_point(data = down.sla1, aes(x = ExpTraitValueSD, y = SpeciesRangeChange), colour = "#0a6933", size = 4, alpha = 0.7) + 
    geom_point(data = up.sla1, aes(x = ExpTraitValueSD, y = SpeciesRangeChange), colour = "#0a6933", size = 7, alpha = 0.7) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    xlab(expression(bold(paste("\nSLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab("75% Relative Species Range Change (%)\n") + range.theme +
    theme(axis.title.x = element_text(face="bold", size=22),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.ticks.length = unit(.25, "cm")))






## HEIGHT VARIATION / RELATIVE CHANGE ----

## Height - 25% quantile model
hei.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitValueSD, data = hei.fut, 
                         iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_25rel_hei_var_mod")
summary(hei.rel25.mod) # positive ns
plot(hei.rel25.mod) 


## Height - median model 
hei.medrel.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, data = hei.fut, 
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_medrel_hei_var_mod")
summary(hei.medrel.mod) # positive ns
plot(hei.medrel.mod) 


## Height - 75% quantile model
hei.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitValueSD, data = hei.fut, 
                      iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_75rel_hei_var_mod")
summary(hei.rel75.mod) # positive ns
plot(hei.rel75.mod) 




## SEED VARIATION / RELATIVE CHANGE ----

## Seed mass - 25% quantile model
seed.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitValueSD, data = seed.fut.all, 
                       iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_25rel_seed_var_mod")
summary(seed.rel25.mod) # positive ns
#plot(seed.rel25.mod) 


## Seed median model
seed.medrel.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, data = seed.fut.all, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_medrel_seed_var_mod")
summary(seed.medrel.mod) # positive ns 
#plot(seed.medrel.mod)


## Seed 75% quantile model
seed.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitValueSD, data = seed.fut.all, 
                       iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_75rel_seed_var_mod")
summary(seed.rel75.mod) # positive ns 
#plot(seed.rel75.mod) 




## FULL MODEL - RELATIVE/VARIATION ----

# Fit full model with all traits 
all.traits.var.rel.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                Hei_var + Seed_var + SLA_var,
                              data = full.traits.yes4, 
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_rel_var_mod")

summary(all.traits.var.rel.mod) # all ns
#plot(all.traits.var.rel.mod)

# Fit full model with all traits and 2x2 interaction 
all.traits.var.rel.mod4 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                 (Hei_var * Seed_var) + (Hei_var * SLA_var) + (SLA_var * Seed_var),
                               data = full.traits.yes4, 
                               iter = 2000, chains = 4, warmup = 400, 
                               file = "models/2022_all_traits_rel_int_var_mod")

summary(all.traits.var.rel.mod4) # all ns
#plot(all.traits.var.rel.mod4)




## SLA VALUES / ABSOLUTE ----

## SLA values - 25% quantile
sla.val.abs25.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitMedianValue, data = sla.fut, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_25abs_sla_val_mod")
summary(sla.val.abs25.mod) #negative ns


## SLA values - median quantile
sla.val.medabs.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = sla.fut, 
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_medabs_sla_val_mod")
summary(sla.val.medabs.mod) #negative ns


# Predictions
slavalabs <- data.frame(LogTraitMedianValue = 
                          seq(from = min(sla.fut$LogTraitMedianValue),
                              to = max(sla.fut$LogTraitMedianValue), by = 0.1))
# predict values
slafit60 <- fitted(
  sla.val.medabs.mod, 
  newdata = slavalabs, 
  re_formula = NULL,
  summary = TRUE
  )

# combine dataframes
colnames(slafit60) = c('fit', 'se', 'lwr', 'upr')
df_sla60 = cbind(slavalabs, slafit60)

# Plot SLA values & ranges 
(sla.absmed.val.plot <- ggplot() + 
    geom_point(data = sla.fut, aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                                   colour = LeafPhenology), size = 5) + 
    scale_colour_manual(values = c("#6CA584","#074923")) +
    xlab(expression(bold(paste("\nSLA values (log mm"^bold("2"), "/mg)")))) +
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    geom_line(data = df_sla60, aes(x = LogTraitMedianValue, y = fit), colour = "#0a6933") + 
    geom_ribbon(data = df_sla60, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), 
                fill = "#0a6933", alpha = 0.5) +
    geom_label_repel(data = subset(sla.fut, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                      "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitMedianValue, cent.abs.median.log, label = SpLabel), color = "black", box.padding = 2, 
                     segment.color = "black", fill = "white", label.size = 1, nudge_x = 0.2,
                     fontface = "italic", size=8) + range.theme + 
    labs(colour = "Deciduousness") +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))


## SLA values - 75% quantile
sla.val.75abs.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitMedianValue, data = sla.fut, 
                          iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_75abs_sla_val_mod")
summary(sla.val.75abs.mod) #negative ns




## SEED VALUES / ABSOLUTE ----

## Seed values - 25% quantile
seed.val.25abs.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitMedianValue, data = seed.fut.all, 
                         iter = 2000, chains = 4, warmup = 400,
                         file = "models/2022_25abs_seed_val_mod")
summary(seed.val.25abs.mod) #negative ns
plot(seed.val.25abs.mod)


# berry species
berry.seed.val.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = berry.spp.rc, 
                             iter = 2000, chains = 4, warmup = 400, 
                             file = "models/2022_seed_val_berry_rc_mod")
summary(berry.seed.val.rc.mod) # negative ns

# wind species
wind.seed.val.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = wind.spp.rc, 
                            iter = 2000, chains = 4, warmup = 400,
                            file = "models/2022_seed_val_wind_rc_mod")
summary(wind.seed.val.rc.mod) # positive ns


# Run the seed mass variation model without the gap-filled species (only true seed mass values)
seed.val.nogap.rc.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = seed.fut, 
                             iter = 2000, chains = 4, warmup = 400,
                             file = "models/2022_seed_val_nogap_rc_mod")
summary(seed.val.nogap.rc.mod) # negative ns




## Seed values - median quantile
seed.val.medabs.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = seed.fut.all, 
                          iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_medabs_seed_val_mod")
summary(seed.val.medabs.mod) # negative ns


# Predictions
seedvalabs <- data.frame(LogTraitMedianValue = 
                          seq(from = min(seed.fut.all$LogTraitMedianValue),
                              to = max(seed.fut.all$LogTraitMedianValue), by = 0.1))
# predict values
seedfit60 <- fitted(
  seed.val.medabs.mod, 
  newdata = seedvalabs, 
  re_formula = NULL,
  summary = TRUE
)

# combine dataframes
colnames(seedfit60) = c('fit', 'se', 'lwr', 'upr')
df_seed60 = cbind(seedvalabs, seedfit60)


# Plotting relationships
(seed.val.medabs.plot <- ggplot() + 
    geom_point(data = other.points.fut,
               aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 5) + 
    geom_point(data = gap.points.fut, shape = 1, 
               aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                   colour = DispersalMode), size = 2.5, stroke = 2) +
    scale_colour_manual(values = c("#CC7000","#FFBA66")) +
    xlab("\nSeed Mass values (log mg)") + 
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    geom_line(data = df_seed60, aes(x = LogTraitMedianValue, y = fit), colour = "#E57E00") + 
    geom_ribbon(data = df_seed60, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr),
                fill = "#E57E00", alpha = 0.4) +
    geom_label_repel(data = subset(seed.fut.all, sp %in% c("Dryas integrifolia", "Linnaea borealis", "Myrica gale")), 
                     aes(LogTraitMedianValue, cent.abs.median.log, label = SpLabel), 
                     color = "black", box.padding = 2, fill = "white",
                     segment.color = "black", label.size = 1,
                     fontface = "italic", size=8) +
    labs(colour = "Dispersal mode") + range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()))


## Seed values - 75% quantile
seed.val.75abs.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitMedianValue, data = seed.fut.all, 
                          iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_75abs_seed_val_mod")
summary(seed.val.75abs.mod) #negative ns
plot(seed.val.75abs.mod)




## HEIGHT VALUES / ABSOLUTE ----

# Plant Height values - 25% quantile
hei.val.abs25.mod <- brm(cent.abs.quan25.log | weights(index) ~ LogTraitMedianValue, data = hei.fut, 
                           iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_25abs_hei_val_mod")
summary(hei.val.abs25.mod) # positive ns



# Plant Height values - median quantile
hei.val.absmed.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, data = hei.fut, 
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_medabs_hei_val_mod")
summary(hei.val.absmed.mod) # positive ns



# Predictions
heivalabs <- data.frame(LogTraitMedianValue = 
                           seq(from = min(hei.fut$LogTraitMedianValue),
                               to = max(hei.fut$LogTraitMedianValue), by = 0.1))
# predict values
heifit60 <- fitted(
  hei.val.absmed.mod, 
  newdata = heivalabs, 
  re_formula = NULL,
  summary = TRUE
)

# combine dataframes
colnames(heifit60) = c('fit', 'se', 'lwr', 'upr')
df_hei60 = cbind(heivalabs, heifit60)


# Plot relationships
(hei.val.medabs.plot <- ggplot() + 
    geom_point(data = hei.fut, aes(x = LogTraitMedianValue, y = cent.abs.median.log, 
                                   colour = gf), size = 5) + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) +
    xlab("\nHeight values (log m)") + 
    ylab(expression(bold(paste("Absolute Species Range Shift (log million km\n"^bold("2\n"), ")\n")))) +
    labs(colour = "Functional group") +
    geom_hline(yintercept = -0.04842876, linetype = "dotted", size = 0.5) +
    geom_line(data = df_hei60, aes(x = LogTraitMedianValue, y = fit), colour = "#800080") + 
    geom_ribbon(data = df_hei60, aes(x = LogTraitMedianValue, ymin = lwr, ymax = upr), fill = "#800080", alpha = 0.3) +
    geom_label_repel(data = subset(height.now, sp %in% c("Rhododendron tomentosum", "Dasiphora fruticosa", "Myrica gale",
                                                         "Linnaea borealis", "Cornus sericea", "Dryas integrifolia")), 
                     aes(LogTraitMedianValue, cent.abs.median.log, label = SpLabel), fill = "white",
                     color = "black", box.padding = 2, segment.color = "black",
                     fontface = "italic", size=8, label.size = 1) + range.theme +
    theme(legend.title = element_text(size = 20, face = "bold"), 
          legend.text=element_text(size = 20),
          legend.position = "top", legend.key = element_blank(),
          legend.background = element_blank()) +
    guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE)))



# Plant Height values - 75% quantile
hei.val.abs75.mod <- brm(cent.abs.quan75.log | weights(index) ~ LogTraitMedianValue, data = hei.fut, 
                         iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_75abs_hei_val_mod")
summary(hei.val.abs75.mod) # positive ns




## ABSOLUTE RANGES - TRAIT VALUES & VARIATION PANEL (FIGURE 5) ----
(abs.all.panel <- ggarrange(hei.val.medabs.plot, sla.absmed.val.plot, seed.val.medabs.plot, 
                            hei.medabs.plot, sla.absmed.plot, seed.medabs.plot,
                            labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
                            nrow = 2, ncol = 3, font.label = list(size = 30)))
ggplot2::ggsave(abs.all.panel, filename = "figures/Figure_5.png", 
                width = 60, height = 60, units = "cm", dpi = 500)


## FULL MODEL VALUES/ABSOLUTE ----

# convert 'all.traits' to long format
all.traits.val.fut <- tidyr::pivot_wider(all.traits4, names_from = TraitShort, 
                                      values_from = c(LogTraitMedianValue, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
all.traits.val.fut.short <- all.traits.val.fut %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(SLA_val = collapse(LogTraitMedianValue_SLA), 
                   SeedMass_val = collapse(LogTraitMedianValue_SeedMass), 
                   Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()


# merge dataframes
full.traits.val.fut <- left_join(all.traits.val.fut.short, more.traits4, by = "sp")
full.traits.val.fut.yes <- distinct(full.traits.val.fut, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.val.abs.mod <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                Hei_val + SeedMass_val + SLA_val, data = full.traits.val.fut.yes, 
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_abs_val_mod")

summary(all.traits.val.abs.mod) # all ns



## Fit full model with all traits and interactions 2x2
all.traits.val.abs.mod4 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                 (Hei_val * SeedMass_val) + (Hei_val * SLA_val) + (SeedMass_val * SLA_val),
                               data = full.traits.val.fut.yes, iter = 2000, chains = 4, warmup = 400, 
                               file = "models/2022_all_traits_abs_int_val_mod") 

summary(all.traits.val.abs.mod4) # all ns
#plot(conditional_effects(all.traits.val.abs.mod4), points = TRUE)



## SLA VALUES / RELATIVE ----

## SLA values - 25% quantile
sla.val.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitMedianValue, 
                       data = sla.fut, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_25rel_sla_val_mod")
summary(sla.val.rel25.mod) # negative ns


## SLA values - median quantile
sla.val.relmed.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                       data = sla.fut, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_medrel_sla_val_mod")
summary(sla.val.relmed.mod) # negative ns
 

## SLA values - 75% quantile
sla.val.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitMedianValue, 
                         data = sla.fut, iter = 2000, chains = 4, warmup = 400,
                         file = "models/2022_75rel_sla_val_mod")
summary(sla.val.rel75.mod) # negative ns




## SEED VALUES / RELATIVE ----

## Seed values - 25% quantile
seed.val.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitMedianValue, 
                         data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_25rel_seed_val_mod")
summary(seed.val.rel25.mod) # negative ns



## Seed values - median quantile
seed.val.relmed.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                          data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_medrel_seed_val_mod")
summary(seed.val.relmed.mod) # negative ns



## Seed values - 75% quantile
seed.val.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitMedianValue, 
                          data = seed.fut.all, iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_75rel_seed_val_mod")
summary(seed.val.rel75.mod) # negative ns 




## HEIGHT VALUES / RELATIVE ----

## Height values - 25% quantile
hei.val.rel25.mod <- brm(cent.rel.quan25.log | weights(index) ~ LogTraitMedianValue, 
                          data = hei.fut, iter = 2000, chains = 4, warmup = 400, 
                         control = list(max_treedepth = 11),
                         file = "models/2022_25rel_hei_val_mod")
summary(hei.val.rel25.mod) # negative ns (slope is 0)


## Height values - median quantile
hei.val.relmed.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                           data = hei.fut, iter = 2000, chains = 4, warmup = 400, 
                          file = "models/2022_medrel_hei_val_mod")
summary(hei.val.relmed.mod) # positive ns (slope is 0)



## Height values - 75% quantile
hei.val.rel75.mod <- brm(cent.rel.quan75.log | weights(index) ~ LogTraitMedianValue, 
                          data = hei.fut, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_75rel_hei_val_mod")
summary(hei.val.rel75.mod) # positive ns




## FULL MODEL VALUES/RELATIVE ----

## Fit full model with all traits - with weights
all.traits.val.rel.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                Hei_val + SeedMass_val + SLA_val,
                              data = full.traits.val.fut.yes, iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_rel_val_mod")
summary(all.traits.val.rel.mod) # all ns


## Fit full model with all traits - with weights and 2x2 interactions 
all.traits.val.rel.mod4 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                 (Hei_val * SeedMass_val) + (SLA_val * Hei_val) + (SLA_val * SeedMass_val),
                               data = full.traits.val.fut.yes, iter = 2000, chains = 4, warmup = 400, 
                               file = "models/2022_all_traits_rel_val_int_mod")

summary(all.traits.val.rel.mod4) # all ns
#plot(conditional_effects(all.traits.val.rel.mod4), points = TRUE) 





## GAINS VS LOSSES - ABSOLUTE ----

# load file with range quantile calculations 
ranges.orig <- read.csv("data/range_quan.csv")

# add gain-loss category on the basis of the median range changes
ranges.gain.loss <- ranges.orig %>% 
  mutate(gainloss.abs = ifelse(abs.median <0, "Loss", "Gain")) %>%
  mutate(gainloss.rel = ifelse(rel.median <0, "Loss", "Gain")) %>% 
  dplyr::select(sp, gainloss.abs, gainloss.rel)
                                               
# merge with trait data
sla.ranges <- merge(sla.fut, ranges.gain.loss, by = "sp")
seed.ranges <- merge(seed.fut.all, ranges.gain.loss, by = "sp")
hei.ranges <- merge(hei.fut, ranges.gain.loss, by = "sp")

# filter SLA 
sla.gain.abs <- sla.ranges %>% filter(gainloss.abs == "Gain")
sla.loss.abs <- sla.ranges %>% filter(gainloss.abs == "Loss")

sla.gain.rel <- sla.ranges %>% filter(gainloss.rel == "Gain")
sla.loss.rel <- sla.ranges %>% filter(gainloss.rel == "Loss")

# filter seed mass
seed.gain.abs <- seed.ranges %>% filter(gainloss.abs == "Gain")
seed.loss.abs <- seed.ranges %>% filter(gainloss.abs == "Loss")

seed.gain.rel <- seed.ranges %>% filter(gainloss.rel == "Gain")
seed.loss.rel <- seed.ranges %>% filter(gainloss.rel == "Loss")

# filter height
hei.gain.abs <- hei.ranges %>% filter(gainloss.abs == "Gain")
hei.loss.abs <- hei.ranges %>% filter(gainloss.abs == "Loss")

hei.gain.rel <- hei.ranges %>% filter(gainloss.rel == "Gain")
hei.loss.rel <- hei.ranges %>% filter(gainloss.rel == "Loss")



## SLA variation - gains
sla.var.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                        data = sla.gain.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_sla_var_gain_mod")
summary(sla.var.gain.mod) # positive ns

## SLA variation - losses
sla.var.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                        data = sla.loss.abs, iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_sla_var_loss_mod")
summary(sla.var.loss.mod) # positive ns



## Seed variation - gains
seed.var.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                        data = seed.gain.abs, iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_seed_var_gain_mod")
print(summary(seed.var.gain.mod), digits = 4) # positive significant


# Model predictions
seed.var.gain.data2 = data.frame(LogTraitValueSD = seq(from = min(seed.gain.abs$LogTraitValueSD),
                                                 to = max(seed.gain.abs$LogTraitValueSD), by = 0.1))

fit.seedvargain = fitted(
  seed.var.gain.mod,
  newdata = seed.var.gain.data2,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.seedvargain) = c('fit', 'se', 'lwr', 'upr')
df_seedgain = cbind(seed.var.gain.data2, fit.seedvargain)

## Back-transforming 

# Mean value used for centering data back in file range.full2
mean(range.full2$abs.median.log) #2.577415

# Mean value of the predicted values for Y (fit)
mean(df_seedgain$fit) # 0.3189906

# Back-centre and back-log-transform
df_seedgain_more <- df_seedgain %>% mutate(y_constant = 2.577415+fit) %>% 
  mutate(y_constant_exp = exp(y_constant)) %>% mutate(back_x = exp(LogTraitValueSD))

# Mean value of the un-centred Y values
mean(df_seedgain_more$y_constant) # 2.896406

# Mean value of the un-centred un-logged Y values
mean(df_seedgain_more$y_constant_exp) # 18.27648

# Back-transformed mean value of the un-centred un-logged Y values
exp(2.896406) # 18.10894

# Min abs + 1 for absolute values = 11.531 +1 = 12.531
df_seedgain_more_more <- df_seedgain_more %>% 
  mutate(y_constant_exp_offset = y_constant_exp - 12.531) %>% 
  mutate(y_cnt_exp_offset_mil = y_constant_exp_offset*1000000)

# Fix the issue with the back-transformation here
seed.gain.abs2 <- seed.gain.abs %>% mutate(ExpTraitValueSD = exp(LogTraitValueSD))

# Separate points per weighting
down.seed2 <- filter(seed.gain.abs2, nobs <20)
up.seed2 <- filter(seed.gain.abs2, nobs >=20)

# Plot this
(back.seed.abs.gain.var.plot <- ggplot() + 
    geom_point(data = df_seedgain_more_more, aes(x = back_x, y = y_cnt_exp_offset_mil), colour = "black", size = 4) +
    geom_point(data = down.seed2, aes(x = ExpTraitValueSD, y = AbsoluteRangeChangeKm), colour = "#E57E00", size = 4, alpha = 0.7) + 
    geom_point(data = up.seed2, aes(x = ExpTraitValueSD, y = AbsoluteRangeChangeKm), colour = "#E57E00", size = 7, alpha = 0.7) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    xlab("\nSeed Mass variation (mg)") +
    ylab(expression(bold(paste("Absolute Species Range Expansions (km\n"^bold("2\n"), ")\n")))) + range.theme +
    theme(axis.title.x = element_text(face="bold", size=22),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.ticks.length = unit(.25, "cm")))



## Seed variation - losses
seed.var.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                        data = seed.loss.abs, iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_seed_var_loss_mod")
summary(seed.var.loss.mod) # positive ns 



## Height variation - gains
hei.var.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                         data = hei.gain.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_hei_var_gain_mod")
summary(hei.var.gain.mod) # positive ns 

## Height variation - losses
hei.var.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitValueSD, 
                         data = hei.loss.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_hei_var_loss_mod")
summary(hei.var.loss.mod) # negative ns 



## All traits variation - gains

# combine trait dataframes
all.traits.gains <- rbind(sla.gain.abs, hei.gain.abs, seed.gain.abs)

# convert to long format
all.traits.gains.long <- pivot_wider(all.traits.gains, names_from = TraitShort, 
                                values_from = c(LogTraitValueSD, index))

# define function 
collapse <- function(x) x[!is.na(x)][1]

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.gains.short <- all.traits.gains.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# select the other relevant traits
more.traits.gains <- all.traits.gains %>% 
  dplyr::select(sp, cent.abs.median.log, cent.rel.median.log, gf, 
                Family, Genus, LeafPhenology, DispersalMode, category.abs, range_cat)

# merge dataframes
full.traits.gains <- left_join(all.traits.gains.short, more.traits.gains, by = "sp")
full.traits.gains.yes <- distinct(full.traits.gains, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.var.abs.gains.mod <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                Hei_var + Seed_var + SLA_var,
                              data = full.traits.gains.yes, iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_abs_var_gain_mod")
summary(all.traits.var.abs.gains.mod) # all ns
#plot(all.traits.var.abs.gains.mod)


## Fit full model with all traits with interactions 2x2
all.traits.var.abs.gains.mod3 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                       (Hei_var * Seed_var) + (Hei_var * SLA_var) + (Seed_var * SLA_var),
                                     data = full.traits.gains.yes, iter = 2000, chains = 4, warmup = 400, 
                                     file = "models/2022_all_traits_abs_var_int_gain_mod")
summary(all.traits.var.abs.gains.mod3) # ns for all
#plot(all.traits.var.abs.gains.mod3)



## All traits variation - losses

# combine trait dataframes
all.traits.losses <- rbind(sla.loss.abs, hei.loss.abs, seed.loss.abs)

# convert to long format
all.traits.losses.long <- pivot_wider(all.traits.losses, names_from = TraitShort, 
                                     values_from = c(LogTraitValueSD, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.losses.short <- all.traits.losses.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# select the other relevant traits
more.traits.losses <- all.traits.losses %>% 
  dplyr::select(sp, cent.abs.median.log, cent.rel.median.log, gf, 
                Family, Genus, LeafPhenology, DispersalMode, category.abs, range_cat)

# merge dataframes
full.traits.losses <- left_join(all.traits.losses.short, more.traits.losses, by = "sp")
full.traits.losses.yes <- distinct(full.traits.losses, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.var.abs.losses.mod <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                      Hei_var + Seed_var + SLA_var,
                                    data = full.traits.losses.yes, iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/2022_all_traits_abs_var_loss_mod")
summary(all.traits.var.abs.losses.mod) # all ns
#plot(all.traits.var.abs.losses.mod)


## Fit full model with all traits and 2x2 interaction
all.traits.var.abs.losses.mod3 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                        (Hei_var * Seed_var) + (SLA_var * Hei_var) + (Seed_var * SLA_var),
                                      data = full.traits.losses.yes, iter = 3000, chains = 4, warmup = 400, 
                                      control = list(adapt_delta = 0.9),
                                      file = "models/2022_all_traits_abs_var_int_loss_mod")
summary(all.traits.var.abs.losses.mod3) # all ns



## SLA values - gains
sla.val.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                        data = sla.gain.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_sla_val_gain_mod")
summary(sla.val.gain.mod) # negative ns

## SLA values - losses
sla.val.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                        data = sla.loss.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_sla_val_loss_mod")
summary(sla.val.loss.mod) # negative ns



## Seed values - gains
seed.val.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                         data = seed.gain.abs, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_seed_val_gain_mod")
summary(seed.val.gain.mod) # negative ns

## Seed values - losses
seed.val.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                         data = seed.loss.abs, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_seed_val_loss_mod")
summary(seed.val.loss.mod) # negative ns



## Height values - gains
hei.val.gain.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                        data = hei.gain.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_hei_val_gain_mod")
summary(hei.val.gain.mod) # negative ns

## Height values - losses
hei.val.loss.mod <- brm(cent.abs.median.log | weights(index) ~ LogTraitMedianValue, 
                        data = hei.loss.abs, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_hei_val_loss_mod")
summary(hei.val.loss.mod) # negative ns



## All traits values - gains

# convert to long format
all.traits.gains.val.long <- pivot_wider(all.traits.gains, names_from = TraitShort, 
                                     values_from = c(LogTraitMedianValue, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.gains.val.short <- all.traits.gains.val.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   Seed_val = collapse(LogTraitMedianValue_SeedMass), 
                   SLA_val = collapse(LogTraitMedianValue_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# merge dataframes
full.traits.val.gains <- left_join(all.traits.gains.val.short, more.traits.gains, by = "sp")
full.traits.gains.val.yes <- distinct(full.traits.val.gains, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.val.abs.gains.mod <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                      Hei_val + Seed_val + SLA_val, data = full.traits.gains.val.yes, 
                                    iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/2022_all_traits_val_abs_gain_mod")
summary(all.traits.val.abs.gains.mod) # all ns


## Fit full model with all traits and interaction
all.traits.val.abs.gains.mod3 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                       (Hei_val * Seed_val) + (Seed_val * SLA_val) + (Hei_val * SLA_val),
                                     data = full.traits.gains.val.yes, 
                                     iter = 2000, chains = 4, warmup = 400, 
                                     file = "models/2022_all_traits_val_abs_int_gain_mod")
summary(all.traits.val.abs.gains.mod3) # all ns



## All traits values - losses

# convert to long format
all.traits.losses.val.long <- pivot_wider(all.traits.losses, names_from = TraitShort, 
                                      values_from = c(LogTraitMedianValue, index))

# merge rows per species and remove species without all 3 trait values (n = 34)
# calculate combined index by normalizing the individual trait indexes
all.traits.losses.val.short <- all.traits.losses.val.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   Seed_val = collapse(LogTraitMedianValue_SeedMass), 
                   SLA_val = collapse(LogTraitMedianValue_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()

# merge dataframes
full.traits.losses.val <- left_join(all.traits.losses.val.short, more.traits.losses, by = "sp")
full.traits.losses.val.yes <- distinct(full.traits.losses.val, sp, .keep_all = TRUE)


## Fit full model with all traits
all.traits.val.abs.losses.mod <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                       Hei_val + Seed_val + SLA_val, data = full.traits.losses.val.yes, 
                                     iter = 2000, chains = 4, warmup = 400, 
                                     file = "models/2022_all_traits_val_abs_loss_mod")
summary(all.traits.val.abs.losses.mod) # seed value negative significant


## Fit full model with all traits with 2x2 interaction 
all.traits.val.abs.losses.mod3 <- brm(cent.abs.median.log | weights(combined.index) ~ 
                                        (Hei_val * Seed_val) + (Hei_val * SLA_val) + (Seed_val * SLA_val),
                                      data = full.traits.losses.val.yes, iter = 2000, chains = 4, warmup = 400, 
                                      control = list(adapt_delta = 0.99, max_treedepth = 11),
                                      file = "models/2022_all_traits_val_abs_int_loss_mod")
summary(all.traits.val.abs.losses.mod3) # ns
plot(conditional_effects(all.traits.val.abs.losses.mod3), points = TRUE) 





## GAINS VS LOSSES - RELATIVE ----

## SLA variation - gains
sla.var.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                        data = sla.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_sla_var_rel_gain_mod")
summary(sla.var.rel.gain.mod) # significant positive 


# Model predictions
sla.var.gain.rel.data = data.frame(LogTraitValueSD = seq(from = min(sla.gain.rel$LogTraitValueSD),
                                                  to = max(sla.gain.rel$LogTraitValueSD), by = 0.05))

fit.slagainrel = fitted(
  sla.var.rel.gain.mod,
  newdata = sla.var.gain.rel.data,
  re_formula = NULL, # ignore random effects
  summary = TRUE # mean and 95% CI
)

colnames(fit.slagainrel) = c('fit', 'se', 'lwr', 'upr')
df_slarelgain = cbind(sla.var.gain.rel.data, fit.slagainrel)


## Back-transforming 

# Mean value used for centering data back in file range.full2
mean(range.full2$rel.median.log) #0.7195768

# Mean value of the predicted values for Y (fit)
mean(df_slarelgain$fit) # 0.2878266

# Back-centre and back-log-transform
df_slarelgain_more <- df_slarelgain %>% mutate(y_constant = 0.7195768+fit) %>% 
  mutate(y_constant_exp = exp(y_constant)) %>% mutate(back_x = exp(LogTraitValueSD))

# Mean value of the un-centred Y values
mean(df_slarelgain_more$y_constant) # 1.007403

# Mean value of the un-centred un-logged Y values
mean(df_slarelgain_more$y_constant_exp) # 2.760432

# Min abs + 1 for relative values = 1+1 = 2
df_slarelgain_more_more <- df_slarelgain_more %>% 
  mutate(y_constant_exp_offset = y_constant_exp - 2) %>% 
  mutate(y_cnt_exp_offset_mil = y_constant_exp_offset*100)

# Fix the issue with the back-transformation here
sla.var.gain.rel2 <- sla.gain.rel %>% mutate(ExpTraitValueSD = exp(LogTraitValueSD))

# Separate points per weighting
down.sla4 <- filter(sla.var.gain.rel2, nobs <20)
up.sla4 <- filter(sla.var.gain.rel2, nobs >=20)


# Plot the right version with the correct back-transform of X data
(back.sla.rel.gain.var.plot <- ggplot() + 
    geom_point(data = df_slarelgain_more_more, aes(x = back_x, y = y_cnt_exp_offset_mil), colour = "black", size = 4) +
    geom_point(data = down.sla4, aes(x = ExpTraitValueSD, y = SpeciesRangeChange), colour = "#0a6933", size = 4, alpha = 0.7) + 
    geom_point(data = up.sla4, aes(x = ExpTraitValueSD, y = SpeciesRangeChange), colour = "#0a6933", size = 7, alpha = 0.7) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    xlab(expression(bold(paste("SLA variation (log mm"^bold("2"), "/mg)")))) +
    ylab("Relative Species Range Expansions (%)\n") + range.theme +
    theme(axis.title.x = element_text(face="bold", size=22),
          axis.text.x  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.title.y = element_text(face="bold", size=22),
          axis.text.y  = element_text(vjust=0.5, size=22, colour = "black"), 
          axis.ticks.length = unit(.25, "cm")))



## SLA variation - losses
sla.var.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                        data = sla.loss.rel, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_sla_var_rel_loss_mod")
summary(sla.var.rel.loss.mod) # positive ns 




## Seed variation - gains
seed.var.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                         data = seed.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_seed_var_rel_gain_mod")
summary(seed.var.rel.gain.mod) # positive significant

## Seed variation - losses
seed.var.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                         data = seed.loss.rel, iter = 2000, chains = 4, warmup = 400, 
                         file = "models/2022_seed_var_rel_loss_mod")
summary(seed.var.rel.loss.mod) # positive ns



## Height variation - gains
hei.var.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                        data = hei.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_hei_var_rel_gain_mod")
summary(hei.var.rel.gain.mod) # positive ns

## Height variation - losses
hei.var.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitValueSD, 
                        data = hei.loss.rel, iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_hei_var_rel_loss_mod")
summary(hei.var.rel.loss.mod) # positive ns




## All traits variation relative - gains

# combine trait dataframes
all.traits.rel.gains <- rbind(sla.gain.rel, hei.gain.rel, seed.gain.rel)

# convert to long format
all.traits.rel.gains.long <- pivot_wider(all.traits.rel.gains, names_from = TraitShort, 
                                     values_from = c(LogTraitValueSD, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.rel.gains.short <- all.traits.rel.gains.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% 
  ungroup()


# merge dataframes
full.traits.rel.gains <- left_join(all.traits.rel.gains.short, more.traits.gains, by = "sp")
full.traits.rel.gains.yes <- distinct(full.traits.rel.gains, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.var.rel.gains.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                      Hei_var + Seed_var + SLA_var, data = full.traits.rel.gains.yes, 
                                    iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/2022_all_traits_var_rel_gain_mod")
summary(all.traits.var.rel.gains.mod) # all ns


## Fit full model with all traits and interaction 2x2 
all.traits.var.rel.gains.mod3 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                       (Hei_var * Seed_var) + (SLA_var * Seed_var) + (Hei_var * SLA_var),
                                     data = full.traits.rel.gains.yes, iter = 2000, chains = 4, warmup = 400, 
                                     control = list(adapt_delta = 0.95), 
                                    file = "models/2022_all_traits_var_int_rel_gain_mod")
summary(all.traits.var.rel.gains.mod3) # all ns
#plot(conditional_effects(all.traits.var.rel.gains.mod3), points = TRUE)




## All traits variation - losses

# combine trait dataframes
all.traits.rel.loss <- rbind(sla.loss.rel, hei.loss.rel, seed.loss.rel)

# convert to long format
all.traits.losses.rel.long <- pivot_wider(all.traits.rel.loss, names_from = TraitShort, 
                                      values_from = c(LogTraitValueSD, index))

# merge rows per species and remove species without all 3 trait values (n = 34)
# calculate combined index by normalizing the individual trait indexes
all.traits.losses.rel.short <- all.traits.losses.rel.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_var = collapse(LogTraitValueSD_PlantHeight), 
                   Seed_var = collapse(LogTraitValueSD_SeedMass), 
                   SLA_var = collapse(LogTraitValueSD_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% ungroup()


# merge dataframes
full.traits.rel.losses <- left_join(all.traits.losses.rel.short, more.traits.losses, by = "sp")
full.traits.rel.losses.yes <- distinct(full.traits.rel.losses, sp, .keep_all = TRUE)


## Fit full model with all traits
all.traits.var.rel.losses.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                       Hei_var + Seed_var + SLA_var, data = full.traits.rel.losses.yes, 
                                     iter = 2000, chains = 4, warmup = 400, 
                                     file = "models/2022_all_traits_var_rel_loss_mod")
summary(all.traits.var.rel.losses.mod) # all ns



## Fit full model with all traits with interaction 2x2 
all.traits.var.rel.losses.mod3 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                        (Hei_var * Seed_var) + (SLA_var * Seed_var) + (Hei_var * SLA_var),
                                      data = full.traits.rel.losses.yes, iter = 2000, chains = 4, warmup = 400, 
                                      control = list(adapt_delta = 0.95), 
                                      file = "models/2022_all_traits_var_int_rel_loss_mod")
summary(all.traits.var.rel.losses.mod3) # all ns






## SLA values - gains
sla.val.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                            data = sla.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/2022_sla_val_rel_gain_mod")
summary(sla.val.rel.gain.mod) # negative ns

## SLA values - losses
sla.val.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                            data = sla.loss.rel, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/2022_sla_val_rel_loss_mod")
summary(sla.val.rel.loss.mod) # positive ns



## Seed values - gains
seed.val.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                             data = seed.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                             file = "models/2022_seed_val_rel_gain_mod")
summary(seed.val.rel.gain.mod) # negative ns

## Seed values - losses
seed.val.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                             data = seed.loss.rel, iter = 2000, chains = 4, warmup = 400, 
                             file = "models/2022_seed_val_rel_loss_mod")
summary(seed.val.rel.loss.mod) # negative ns



## Height values - gains
hei.val.rel.gain.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                            data = hei.gain.rel, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/2022_hei_val_rel_gain_mod")
summary(hei.val.rel.gain.mod) # negative ns

## Height values - losses
hei.val.rel.loss.mod <- brm(cent.rel.median.log | weights(index) ~ LogTraitMedianValue, 
                            data = hei.loss.rel, iter = 2000, chains = 4, warmup = 400, 
                            file = "models/2022_hei_val_rel_loss_mod")
summary(hei.val.rel.loss.mod) # negative ns



## All traits values relative - gains

# convert to long format
all.traits.rel.val.gains.long <- pivot_wider(all.traits.rel.gains, names_from = TraitShort, 
                                         values_from = c(LogTraitMedianValue, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.rel.val.gains.short <- all.traits.rel.val.gains.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   Seed_val = collapse(LogTraitMedianValue_SeedMass), 
                   SLA_val = collapse(LogTraitMedianValue_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% 
  ungroup()


# merge dataframes
full.traits.rel.val.gains <- left_join(all.traits.rel.val.gains.short, more.traits.gains, by = "sp")
full.traits.rel.val.gains.yes <- distinct(full.traits.rel.val.gains, sp, .keep_all = TRUE)


## Fit full model with all traits 
all.traits.val.rel.gains.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                      Hei_val + Seed_val + SLA_val, data = full.traits.rel.val.gains.yes, 
                                    iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/2022_all_traits_val_rel_gain_mod")
summary(all.traits.val.rel.gains.mod) # all ns



## Fit full model with all traits with interaction 2x2 
all.traits.val.rel.gains.mod3 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                       (Hei_val * Seed_val) + (SLA_val * Seed_val) + (Hei_val * SLA_val),
                                     data = full.traits.rel.val.gains.yes, iter = 2000, chains = 4, warmup = 400, 
                                     control = list(adapt_delta = 0.9),
                                     file = "models/2022_all_traits_val_int_rel_gain_mod")
summary(all.traits.val.rel.gains.mod3) # all ns




## All traits values relative - losses

# convert to long format
all.traits.rel.val.loss.long <- pivot_wider(all.traits.rel.loss, names_from = TraitShort, 
                                             values_from = c(LogTraitMedianValue, index))

# merge rows per species and remove species without all 3 trait values (n = 35)
# calculate combined index by normalizing the individual trait indexes
all.traits.rel.val.loss.short <- all.traits.rel.val.loss.long %>% dplyr::group_by(sp) %>% 
  dplyr::summarise(Hei_val = collapse(LogTraitMedianValue_PlantHeight), 
                   Seed_val = collapse(LogTraitMedianValue_SeedMass), 
                   SLA_val = collapse(LogTraitMedianValue_SLA),
                   SLA_Index = collapse(index_SLA), 
                   SeedMass_Index = collapse(index_SeedMass),
                   Hei_Index = collapse(index_PlantHeight)) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(combined.index = (SLA_Index + SeedMass_Index + Hei_Index)/3) %>% 
  ungroup()


# merge dataframes
full.traits.rel.val.loss <- left_join(all.traits.rel.val.loss.short, more.traits.losses, by = "sp")
full.traits.rel.val.loss.yes <- distinct(full.traits.rel.val.loss, sp, .keep_all = TRUE)


## Fit full model with all traits - with weights
all.traits.val.rel.loss.mod <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                      Hei_val + Seed_val + SLA_val, data = full.traits.rel.val.loss.yes, 
                                    iter = 2000, chains = 4, warmup = 400, 
                                   file = "models/2022_all_traits_val_rel_loss_mod")
summary(all.traits.val.rel.loss.mod) # all ns



## Fit full model with all traits - with weights & interaction 2x2
all.traits.val.rel.loss.mod3 <- brm(cent.rel.median.log | weights(combined.index) ~ 
                                      (Hei_val * Seed_val) + (SLA_val * Seed_val) + (Hei_val * SLA_val),
                                    data = full.traits.rel.val.loss.yes, iter = 2000, chains = 4, warmup = 400, 
                                    control = list(adapt_delta = 0.97, max_treedepth = 12), 
                                    file = "models/2022_all_traits_val_int_rel_loss_mod")
summary(all.traits.val.rel.loss.mod3) # SLA positive significant, Height x SLA positive significant
plot(conditional_effects(all.traits.val.rel.loss.mod3), points = TRUE)





## CATEGORICAL MODELS ----

## Absolute range changes & functional group
future.fg.mod <- brm(cent.abs.median.log ~ gf, data = full.traits.yes4, 
                      iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_abs_fg_mod")
summary(future.fg.mod) # no significant difference between the three functional groups
conditional_effects(future.fg.mod)


## Absolute range changes & dispersal mode
future.disp.mod <- brm(cent.abs.median.log ~ DispersalMode, data = full.traits.yes4, 
                        iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_abs_disp_mod")
summary(future.disp.mod) # no significant difference between the two dispersal modes
conditional_effects(future.disp.mod)


## Absolute range changes & decidiousness
future.decid.mod <- brm(cent.abs.median.log ~ LeafPhenology, data = full.traits.yes4, 
                         iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_abs_decid_mod")
summary(future.decid.mod) # no significant difference between evergreen and deciduous
conditional_effects(future.decid.mod)


## Absolute range changes & family
future.fam.mod <- brm(cent.abs.median.log ~ Family, data = full.traits.yes4, 
                       iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_abs_fam_mod")
summary(future.fam.mod) # no significant difference between the 7 families
conditional_effects(future.fam.mod)



## Relative range changes & functional group
rel.fg.mod <- brm(cent.rel.median.log ~ gf, data = full.traits.yes4, 
                     iter = 2000, chains = 4, warmup = 400,
                  file = "models/2022_rel_fg_mod")
summary(rel.fg.mod) # no significant difference between the three functional groups
conditional_effects(rel.fg.mod)


## Relative range changes & dispersal mode
rel.disp.mod <- brm(cent.rel.median.log ~ DispersalMode, data = full.traits.yes4, 
                       iter = 2000, chains = 4, warmup = 400,
                    file = "models/2022_rel_disp_mod")
summary(rel.disp.mod) # no significant difference between the two dispersal modes
conditional_effects(rel.disp.mod)


## Relative range changes & decidiousness
rel.decid.mod <- brm(cent.rel.median.log ~ LeafPhenology, data = full.traits.yes4, 
                        iter = 2000, chains = 4, warmup = 400, 
                     file = "models/2022_rel_decid_mod")
summary(rel.decid.mod) # no significant difference between evergreen and deciduous
conditional_effects(rel.decid.mod)


## Relative range changes & family
rel.fam.mod <- brm(cent.rel.median.log ~ Family, data = full.traits.yes4, 
                      iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_rel_fam_mod")
summary(rel.fam.mod) # no significant difference between the 7 families
conditional_effects(rel.fam.mod)



## PANEL WITH BACK-TRANSFORMATIONS ----
(backtrans.panel <- ggarrange(back.sla.abs75.var.plot2,
                              back.sla.rel75.var.plot2, back.sla.rel.gain.var.plot,
                              back.seed.abs.gain.var.plot,
                            nrow = 2, ncol = 2, 
                            labels = c("(a)", "(b)", "(c)", "(d)"),
                            font.label = list(size = 30)))

ggplot2::ggsave(backtrans.panel, filename = "scripts/users/mgarciacriado/traits_vs_ranges/figures/2022_FigS4_Panel.png", 
                width = 50, height = 50, units = "cm")






## COVER CHANGE MODELS ----
cover.slopes <- read.csv("data/CoverChangeSpecies_QuantsCats_Final.csv")
names(cover.slopes)[names(cover.slopes) == 'Name'] <- 'sp' # (36 spps with cover data)


# join the dataframes with trait & range info with the cover change slopes
seed.cov <- left_join(seed.fut.all, cover.slopes, by = "sp")
sla.cov <- left_join(sla.fut, cover.slopes, by = "sp")
hei.cov <- left_join(hei.fut, cover.slopes, by = "sp")


#### SLA VARIATION MODEL
sla.cov.var.mod <- brm(MeanSlope | weights(index) ~ LogTraitValueSD, data = sla.cov,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/2022_sla_var_cov_mod")
summary(sla.cov.var.mod) # negative ns



#### SEED VARIATION MODEL
seed.cov.var.mod <- brm(MeanSlope | weights(index) ~ LogTraitValueSD, data = seed.cov,
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_seed_var_cov_mod")
summary(seed.cov.var.mod) # positive ns



#### HEIGHT VARIATION MODEL
hei.cov.var.mod <- brm(MeanSlope | weights(index) ~ LogTraitValueSD, data = hei.cov,
                        iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_hei_var_cov_mod")
summary(hei.cov.var.mod) # positive ns



#### ALL VARIATION MODEL
all.traits.var.cov.df <- left_join(full.traits.yes4, cover.slopes, by = "sp")
all.traits.var.cov.df2 <- all.traits.var.cov.df %>% drop_na(MeanSlope)

## Fit full model with all traits variation
all.traits.var.cov.mod <- brm(MeanSlope | weights(combined.index) ~ 
                                Hei_var + Seed_var + SLA_var, data = all.traits.var.cov.df2,
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_var_cov_mod")
summary(all.traits.var.cov.mod) # ns


#### SLA VALUES MODEL
sla.cov.val.mod <- brm(MeanSlope | weights(index) ~ LogTraitMedianValue, data = hei.cov, 
                          iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_sla_val_cov_mod")

summary(sla.cov.val.mod) # positive ns



#### SEED VALUES MODEL
seed.cov.val.mod <- brm(MeanSlope | weights(index) ~ LogTraitMedianValue, data = seed.cov, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_seed_val_cov_mod")

summary(seed.cov.val.mod) # negative ns



#### HEIGHT VALUES MODEL
hei.cov.val.mod <- brm(MeanSlope | weights(index) ~ LogTraitMedianValue, data = hei.cov, 
                        iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_hei_val_cov_mod")

summary(hei.cov.val.mod) # positive ns



#### ALL TRAITS VALUES MODEL
all.traits.val.cov.df <- left_join(full.traits.val.fut.yes, cover.slopes, by = "sp")
all.traits.val.cov.df2 <- all.traits.val.cov.df %>% drop_na(MeanSlope)

## Fit full model with all traits variation - with weights
all.traits.val.cov.mod <- brm(MeanSlope | weights(combined.index) ~ 
                                Hei_val + SeedMass_val + SLA_val,
                              data = all.traits.val.cov.df2,
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/2022_all_traits_val_cov_mod")

summary(all.traits.val.cov.mod) # ns
