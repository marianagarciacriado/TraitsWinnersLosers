## CALCULATE CHANGE IN COVER BY SPECIES PER ITEX SITE ##
## Written by Anne, Stan model by Anne & Nadja, further edits by Mariana ##
## December 2020 ##

# The Beta-Bernoulli models did not work well (too many plots with only 0's for the Bernoulli part?) 
# so instead I went with a Poisson model where cover is between 0 and 100 (rounded up using "ceiling" 
# so any amount <1 becomes 1). This is based on relative cover (not including bryophytes and lichens)



## PACKAGES ----
library(ggplot2)
library(tidyverse)
library(rjags)
library(R2jags)


## DATA IMPORT ----
load("scripts/users/abjorkman/Bayesian_CWM/coverc_sub.RData")

# Exclude morphospecies and some sites and subsites that have only one plot
cover.spp <- coverc.sub[!is.na(coverc.sub$COVER) & coverc.sub$Morphosp==0 & 
                          coverc.sub$SITE %in% c("TAISETSU","TIBET","STEPSTONES","FAROE","SADVENT")==F & 
                          coverc.sub$SUBSITE %in% c("ALEXFIORD:CASSIOPE_COVER","BARROW:ARCTOPHILA_POND_ORDINATION")==F & 
                          coverc.sub$TRTMT == "CTL", ]

# Because Iceland has so many sites, temporarily make each subsite within a site the same name (to calculate num of sites)
cover.spp$SUBSITEtemp <- cover.spp$SUBSITE
cover.spp$SUBSITEtemp[cover.spp$SITE %in% c("AKUREYRI","AUDKULUHEIDI","BLONDUOS","DALSMYNNI","HJARDARLAND",
                                            "HOLTAVORDUHEIDI","MODRUVELLIR","OXNADALSHEIDI","THINGVELLIR",
                                            "THYKKVIBAER")] <- 
  cover.spp$SITE[cover.spp$SITE %in% c("AKUREYRI","AUDKULUHEIDI","BLONDUOS","DALSMYNNI","HJARDARLAND",
                                       "HOLTAVORDUHEIDI","MODRUVELLIR","OXNADALSHEIDI","THINGVELLIR","THYKKVIBAER")]


# SUBSET TO MARIANA'S SPECIES 
mar.spp <- c("Dasiphora fruticosa", "Myrica gale", "Rhododendron lapponicum", "Ribes nigrum", "Ribes spicatum", 
             "Sibbaldia procumbens", "Harrimanella hypnoides", "Arctous rubra", "Salix phylicifolia", 
             "Salix hastata", "Rosa majalis", "Rosa acicularis", "Rubus chamaemorus", "Rubus idaeus", 
             "Artemisia dracunculus", "Salix lanata", "Salix niphoclada", "Artemisia gmelinii", 
             "Salix barrattiana", "Andromeda glaucophylla", "Salix myrsinites", "Salix lapponum", 
             "Vaccinium myrtilloides", "Andromeda polifolia", "Cornus canadensis", "Thymus praecox", 
             "Salix vestita", "Salix argyrocarpa", "Vaccinium oxycoccos", "Betula glandulosa", 
             "Vaccinium caespitosum", "Thymus serpyllum", "Diapensia lapponica", "Betula humilis", 
             "Artemisia frigida", "Ribes glandulosum", "Chamaedaphne calyculata", "Salix rosmarinifolia", 
             "Arctostaphylos uva-ursi", "Cornus sericea", "Linnaea borealis", "Ledum palustre", 
             "Phyllodoce caerulea", "Cassiope tetragona", "Arctous alpina", "Vaccinium myrtillus", 
             "Empetrum nigrum", "Vaccinium vitis-idaea", "Calluna vulgaris", "Salix phlebophylla", 
             "Dryas octopetala", "Salix glauca", "Salix herbacea", "Betula nana", "Salix rotundifolia", 
             "Vaccinium uliginosum", "Salix reticulata", "Loiseleuria procumbens", "Salix arctica", 
             "Salix polaris", "Salix pulchra", "Dryas integrifolia")

mar.spp[mar.spp %in% unique(cover.spp$Name)==F]

# Extra sp: “Phyllodoce caerulea”, “Cassiope tetragona”, “Arctous alpina”, “Vaccinium myrtillus”, “Empetrum nigrum”, “Vaccinium vitis-idaea”, “Calluna vulgaris”, “Salix phlebophylla”, “Dryas octopetala”, “Salix glauca”, “Salix herbacea”, “Betula nana”, “Salix rotundifolia”, “Vaccinium uliginosum”, “Salix reticulata”, “Loiseleuria procumbens”, “Salix Arctica”, “Salix polaris”, “Salix pulchra”, “Dryas integrifolia”

# [1] "Myrica gale"             "Ribes nigrum"            "Ribes spicatum"          "Rosa majalis"           
# [5] "Rosa acicularis"         "Rubus idaeus"            "Artemisia dracunculus"   "Artemisia gmelinii"     
# [9] "Salix barrattiana"       "Andromeda glaucophylla"  "Salix myrsinites"        "Salix lapponum"         
# [13] "Vaccinium myrtilloides"  "Cornus canadensis"       "Salix vestita"           "Salix argyrocarpa"      
# [17] "Vaccinium oxycoccos"     "Vaccinium caespitosum"   "Thymus serpyllum"        "Betula humilis"         
# [21] "Artemisia frigida"       "Ribes glandulosum"       "Chamaedaphne calyculata" "Salix rosmarinifolia"   
# [25] "Cornus sericea"          "Linnaea borealis"  


cover.spp2 <- cover.spp[cover.spp$Name %in% c(mar.spp),] # ALL OF MARIANA'S SPECIES

cover.spp2$SemiUniquePLOT <- paste(cover.spp2$CRUGrid, cover.spp2$SUBSITE, cover.spp2$PLOT, sep="_")

cover.backup <- cover.spp2



# PREPARE FOR MODEL INPUT ----

# Not centering within each Grid Cell because not necessary when slopes estimated separately for each site-species

# For beta distribution, divide by 100 (must be between 0 and 1)
cover.spp2$COVER2 <- cover.spp2$COVER/100

# Species and site numbers MUST be the same in both discrete and continuous datasets 
# so we create the factor before splitting (only relevant for beta-bern)

cover.spp2 <- cover.spp2[order(cover.spp2$Name,cover.spp2$CRUGrid,cover.spp2$SUBSITE, cover.spp2$YEAR, cover.spp2$SemiUniquePLOT),]

cover.spp2$Speciesf <- factor(cover.spp2$Name, levels = unique(cover.spp2$Name))
cover.spp2$Sitef <- factor(cover.spp2$SUBSITE, levels = unique(cover.spp2$SUBSITE)) # USING SUBSITE INSTEAD OF CRU GRID!

cover.spp2$SiteSpeciesName <- paste(cover.spp2$Name,cover.spp2$SUBSITE,sep="_")
cover.spp2$SiteSpeciesYrName <- paste(cover.spp2$SiteSpeciesName, cover.spp2$YEAR, sep="_")

cover.spp2$SiteSpeciesf <- factor(cover.spp2$SiteSpeciesName, levels = unique(cover.spp2$SiteSpeciesName))

cover.spp2$SpeciesNum <- as.numeric(cover.spp2$Speciesf)
cover.spp2$SiteNum <- as.numeric(cover.spp2$Sitef)
cover.spp2$SiteSpecies <- as.numeric(cover.spp2$SiteSpeciesf)

cover.spp2$SiteSppYr<-as.numeric(factor(cover.spp2$SiteSpeciesYrName, levels = unique(cover.spp2$SiteSpeciesYrName)))

# Plot effect = 0 if plots are not repeated OR if the species occurs in only 1 plot in a site
cover.spp2 <- plyr::ddply(cover.spp2, c("Name","SemiUniquePLOT"), transform, # do plots repeat?
               countPlotReps = length(unique(YEAR)))
cover.spp2 <- plyr::ddply(cover.spp2, c("Name","Sitef"), transform, # do plots repeat?
                    countPlotReps2 = length(unique(SemiUniquePLOT)))
cover.spp2$PlotEffect <- ifelse(cover.spp2$countPlotReps > 1, 1, 0)
cover.spp2$PlotEffect <- ifelse(cover.spp2$countPlotReps2 < 2, 0, cover.spp2$PlotEffect)

cover.spp2$SubsiteSpeciesName <- paste(cover.spp2$Name, cover.spp2$SUBSITE, sep="_")
cover.spp2$PlotSpeciesName <- paste(cover.spp2$Name, cover.spp2$SemiUniquePLOT, sep="_")
cover.spp2$PlotSpeciesName[cover.spp2$PlotEffect == 0] <- paste(cover.spp2$SubsiteSpeciesName[cover.spp2$PlotEffect == 0],"plot",sep="_") # give plots that aren't repeating the subsite name since they'll be set to 0 later anyway (so fewer plot effects to estimate - hopefully speed things up)

cover.spp2$PlotSpecies <- as.numeric(factor(cover.spp2$PlotSpeciesName, levels = unique(cover.spp2$PlotSpeciesName)))

# Creating a key with the site/subsite info for later
itex.key <- cover.spp2 %>% select(SUBSITE, SITE) %>% distinct()


# Check raw data
ggplot(cover.spp2)+
  geom_point(aes(x=YEAR, y=COVER2, colour = factor(SiteNum)))+
  geom_smooth(aes(x=YEAR, y=COVER2, colour = factor(SiteNum)), se=F, method = "lm")+
  facet_wrap(~Name, scales = "free_y")+
  theme_bw()

ggplot(cover.spp2)+
  geom_histogram(aes(x=COVER2, fill = factor(SiteNum)), bins = 100)+
  facet_wrap(~Name, scales = "free_y")+
  theme_bw()



#JAGS - Poisson - By PLOT -----


jags.dat<-list(
  Nsppplot = length(unique(cover.spp2$PlotSpecies)),
  #
  Nobs = length(cover.spp2[,1]),
  sppplot = cover.spp2$PlotSpecies,
  covobs = round(ceiling(cover.spp2$COVER2*100),0),
  year = as.numeric(cover.spp2$YEAR - 2003)
)

str(jags.dat)


write("

model {
  
# Likelihood
  
  for (i in 1:Nobs){
    covobs[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha[sppplot[i]] + beta[sppplot[i]]*year[i] + eps[i]
    eps[i]~dnorm(0,tau)
  }
  
# Assess fit
  
for (i in 1:Nobs){
  Presi[i] <- (covobs[i] - lambda[i])/sqrt(lambda[i])
  cov.new[i] ~ dpois(lambda[i])
  Presi.new[i] <- (cov.new[i] - lambda[i])/sqrt(lambda[i])
  D[i] <- pow(Presi[i],2)
  D.new[i] <- pow(Presi.new[i],2)
}

  fit <- sum(D[])
  fit.new <- sum(D.new[])

# Priors

for (i in 1:Nsppplot){
 alpha[i] ~ dnorm(0,0.25) #tau equivalent to a sigma of 2
 beta[i] ~ dnorm(0,0.01)
}

tau <- 1/(sigma*sigma)
sigma ~ dunif(0,10)
  
}


","scripts/users/abjorkman/Cover_Change/cover_change_betabern_poisson.jags")

params <- c("alpha", "beta", "lambda", "Presi", "fit", "fit.new", "sigma")



mout.jags <- jags(jags.dat, inits = NULL, params, 
                  model.file = "scripts/users/abjorkman/Cover_Change/cover_change_betabern_poisson.jags", 
                  n.chains = 3, n.iter = 20000, n.burnin = 15000, n.thin = 2, DIC = FALSE, 
                  working.directory = NULL, progress.bar = "text") 


plot(mout.jags) #check convergence, etc.

# extract coefficients 
coeff.jags<-as.data.frame(mout.jags$BUGSoutput$summary[,c('mean','sd','2.5%','97.5%','Rhat')])

# add identifying info to data frame
coeff.jags$Param <- as.vector(sapply(strsplit(rownames(coeff.jags),"[[]",fixed=FALSE), "[", 1))
coeff.jags$Number <- as.numeric(as.vector(sapply(strsplit(rownames(coeff.jags),"[^0-9]+",fixed=FALSE), "[", 2)))
coeff.jags[coeff.jags$Rhat > 1.1 & !is.na(coeff.jags$Rhat),]

# check fit
plot(mout.jags$BUGSoutput$mean$Presi, las =1)
plot(mout.jags$BUGSoutput$sims.list$fit, mout.jags$BUGSoutput$sims.list$fit.new)

# compare estimates
alphas <- coeff.jags[coeff.jags$Param=="alpha",]
alphas$Plot <- cover.spp2$SemiUniquePLOT[match(alphas$Number,cover.spp2$PlotSpecies)]
alphas$Name <- cover.spp2$Name[match(alphas$Number,cover.spp2$PlotSpecies)]
alphas$SUBSITE <- cover.spp2$SUBSITE[match(alphas$Number,cover.spp2$PlotSpecies)]
alphas$CRUGrid <- cover.spp2$CRUGrid[match(alphas$Number,cover.spp2$PlotSpecies)]
alphas$RepeatedPlots <- cover.spp2$RepeatedPlots[match(alphas$Number,cover.spp2$PlotSpecies)]

betas <- coeff.jags[coeff.jags$Param=="beta",]
betas$Plot <- cover.spp2$SemiUniquePLOT[match(betas$Number,cover.spp2$PlotSpecies)]
betas$Name <- cover.spp2$Name[match(betas$Number,cover.spp2$PlotSpecies)]
betas$SUBSITE <- cover.spp2$SUBSITE[match(betas$Number,cover.spp2$PlotSpecies)]
betas$CRUGrid <- cover.spp2$CRUGrid[match(betas$Number,cover.spp2$PlotSpecies)]
betas$RepeatedPlots <- cover.spp2$RepeatedPlots[match(alphas$Number,cover.spp2$PlotSpecies)]


m <- function(x) {lm(log(round(ceiling(COVER2*100),0)+1) ~ I(YEAR-2003), data = x)}

#safe.m <- dplyr::failwith(NA, m, quiet = TRUE) # if dplyr < 1.0, do this
safe.m <- purrr::possibly(m, otherwise = NA, quiet = TRUE) # if dplyr > 1.0, do this 

mods <- plyr::dlply(cover.spp2, c("PlotSpecies"), safe.m)

outputf.lm <- function (x) c(
  slopes = summary(x)$coef[2,1],
  intercepts = summary(x)$coef[1,1]
)

#safe.out <- dplyr::failwith(NA, outputf.lm, quiet = T) # if dplyr < 1.0, do this
safe.out <- purrr::possibly(outputf.lm, otherwise = NA, quiet = TRUE) # if dplyr > 1.0, do this 

mods.out <- plyr::ldply(mods,safe.out)

# compare modeled and raw data
compare <- merge(betas, mods.out, by.x="Number", by.y ="PlotSpecies")
compare.alphas <- merge(alphas, mods.out, by.x="Number", by.y ="PlotSpecies")

ggplot(compare[!is.na(compare$slopes),])+
  geom_hline(yintercept=0, linetype = "dotted")+
  geom_vline(xintercept=0, linetype = "dotted")+
  geom_point(aes(x=mean,y=slopes)) 

ggplot(compare.alphas[!is.na(compare.alphas$intercepts),])+
  geom_hline(yintercept=0, linetype = "dotted")+
  geom_vline(xintercept=0, linetype = "dotted")+
  geom_point(aes(x=exp(mean),y=exp(intercepts)))

ggplot(cover.spp2[cover.spp2$SemiUniquePLOT=="71429_NIWOT:MOIST MEADOW_SADDLE_31" & cover.spp2$Name=="Sibbaldia procumbens",])+
  geom_point(aes(x=YEAR, y = COVER2*100))+
  geom_smooth(aes(x=YEAR, y = COVER2*100))

ggplot(cover.spp2[cover.spp2$SUBSITE=="NIWOT:MOIST MEADOW_SADDLE" & cover.spp2$Name=="Sibbaldia procumbens",])+
  geom_point(aes(x=YEAR-2003, y = log(COVER2*100+1), colour = SemiUniquePLOT))+
  geom_smooth(aes(x=YEAR-2003, y = log(COVER2*100+1), colour = SemiUniquePLOT), se = F, method="lm")


# JAGS output - derived parameters ----

# Extract iterations

iter.alpha <- as.data.frame(t(mout.jags$BUGSoutput$sims.list$alpha), stringsAsFactors = F)
iter.alpha$PlotSpecies <- 1:length(iter.alpha[,1])
iter.alpha$Name <- cover.spp2$Name[match(iter.alpha$PlotSpecies,cover.spp2$PlotSpecies)]
iter.alpha$SUBSITE <- cover.spp2$SUBSITE[match(iter.alpha$PlotSpecies,cover.spp2$PlotSpecies)]

# If dplyr < 1.0, do this: convert to long format and add Site info
iter.alpha.long <- iter.alpha %>% pivot_longer(cols = starts_with("V"), names_to = "iterations", values_to = "alpha")
iter.site <- merge(iter.alpha.long, itex.key, by = "SUBSITE")

iter.beta <- as.data.frame(t(mout.jags$BUGSoutput$sims.list$beta), stringsAsFactors = F)
iter.beta$PlotSpecies <- 1:length(iter.beta[,1])
iter.beta$Name <- cover.spp2$Name[match(iter.beta$PlotSpecies,cover.spp2$PlotSpecies)]
iter.beta$SUBSITE <- cover.spp2$SUBSITE[match(iter.beta$PlotSpecies,cover.spp2$PlotSpecies)]

save(iter.beta, file = "scripts/users/mgarciacriado/traits_vs_ranges/data/iter_beta.RData")

# If dplyr < 1.0, do this: convert to long format and add Site info
iter.beta.long <- iter.beta %>% pivot_longer(cols = starts_with("V"), names_to = "iterations", values_to = "beta")
iter.site.beta <- merge(iter.beta.long, itex.key, by = "SUBSITE")




# Mean of all iterations within a SUBSITE (will be the same for non-repeating plot subsites)
mean.alpha <- iter.alpha %>%
  dplyr::group_by(Name, SUBSITE) %>%
  dplyr::summarise(dplyr::across(everything(), list(mean)))


## We tried to do this a different way since across didn't work on my computer
## Aggregating so we don't lose subsite and site variability
mean.alpha.subsite <- iter.site %>% 
  group_by(Name, SUBSITE) %>% mutate(MeanSubsite = median(alpha)) %>% select(Name, SITE, SUBSITE, MeanSubsite) %>% distinct()

mean.alpha.site <- mean.alpha.subsite %>% 
  group_by(Name, SITE) %>% summarise(MeanSite = median(MeanSubsite)) 

mean.alpha.sp <- mean.alpha.site %>% 
  group_by(Name) %>% summarise(MeanSpecies = median(MeanSite)) 



## betas
mean.beta <- iter.beta %>%
  group_by( Name, SUBSITE ) %>%
  dplyr::summarise(across(everything(), list(mean)))


# alternative way as above
mean.beta.subsite <- iter.site.beta %>% 
  group_by(Name, SUBSITE) %>% mutate(MeanSubsite = mean(beta)) %>% mutate(MedianSubsite = median(beta)) %>%
  select(Name, SITE, SUBSITE, MeanSubsite, MedianSubsite) %>% distinct()

mean.beta.site <- mean.beta.subsite %>% 
  group_by(Name, SITE) %>% summarise(MeanSite = mean(MeanSubsite), MedianSite = median(MedianSubsite))

mean.beta.sp <- mean.beta.site %>% 
  group_by(Name) %>% summarise(MeanSpecies = mean(MeanSite), MedianSpecies = median(MedianSite))



# Anne's code continues here
mean.alpha[1:10,1:10]
mean.beta[1:10,1:10]

# Quantiles per subsite
quant.alpha <- bind_cols(as.data.frame(t(apply(as.matrix(mean.alpha[,3:ncol(mean.alpha)]), 1, quantile, c(0.05, 0.5, 0.95)))))
quant.alpha$Name <- mean.alpha$Name
quant.alpha$SUBSITE <- mean.alpha$SUBSITE

quant.beta <- bind_cols(as.data.frame(t(apply(as.matrix(mean.beta[,3:ncol(mean.beta)]), 1, quantile, c(0.05, 0.5, 0.95)))))
quant.beta$Name <- mean.beta$Name
quant.beta$SUBSITE <- mean.beta$SUBSITE

# Compare to raw data (spot check)
ggplot(cover.spp2[cover.spp2$Name=="Arctostaphylos uva-ursi",])+
  geom_point(aes(x=YEAR-2003, y = log(COVER2*100+1)), position = "jitter")+
  geom_smooth(aes(x=YEAR-2003, y = log(COVER2*100+1)), se = F, method="lm")+
  facet_wrap(~SUBSITE, scales = "free")

# Save output

write.csv(quant.beta, file = "CoverChange_MarianaSpp_bySubsite.csv")
write.csv(quant.beta, file = "scripts/users/mgarciacriado/traits_vs_ranges/data/CoverChange_AllSpecies.csv")



# Let's do the same as above for all species
species.beta <- mean.beta %>%
  group_by(Name) %>%
  dplyr::summarise(across(everything(), list(mean)))

# Calculate quantiles and put this in a consistent format
all.beta <- bind_cols(as.data.frame(t(apply(as.matrix(species.beta[,3:ncol(species.beta)]), 1, quantile, c(0.05, 0.5, 0.95)))))
all.beta$Name <- species.beta$Name
colnames(all.beta) <- c("LowerQuantile", "MeanSlope", "UpperQuantile", "Name")

# Calculate species categories based on cover change
all.cover.quan.cat <- all.beta %>% mutate(CoverCategory = case_when(LowerQuantile > 0 ~ "Winner",
                                   UpperQuantile < 0 ~ "Loser",
                                   TRUE ~ "No change"))

# Save this file
write.csv(all.cover.quan.cat, "data/CoverChangeSpecies_QuantsCats_Final.csv")









