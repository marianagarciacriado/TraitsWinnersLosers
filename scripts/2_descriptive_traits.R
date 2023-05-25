## Trait-range manuscript
## Mariana Garcia
## February 2020
## Script 2. Descriptive trait data


## LIBRARIES ----
library(tidyverse)
library(ggpubr)
library(ggOceanMaps)
library(cowplot)

#devtools::install_github("MikkoVihtakari/ggOceanMaps")
#devtools::install_github("MikkoVihtakari/ggOceanMapsData")


## THEME ----
range.theme <- theme(legend.position = "none",
                     axis.title.x = element_text(face="bold", size=19.5),
                     axis.text.x  = element_text(vjust=0.5, size=16, colour = "black"), 
                     axis.title.y = element_text(face="bold", size=19.5),
                     axis.text.y  = element_text(vjust=0.5, size=16, colour = "black"),
                     panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
                     panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                     plot.title = element_text(color = "black", size = 18, face = "bold", hjust = 0.5),
                     plot.margin = unit(c(1,1,1,1), units = , "cm"))


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## LOAD & FILTER DATA ----
load("data/2022_master_new_superfinal.RData")
#this is assigned to the object 'master.new.superfinal' - 17921 obs




## DATAFRAMES PER TRAIT ----
sla.bubble <- master.new.superfinal %>% filter(TraitShort == "SLA") %>% 
  group_by(sp) %>% mutate(MedianValue = median(StdValue)) %>% ungroup()

seed.bubble <- master.new.superfinal %>% filter(TraitShort == "SeedMass") %>%
  group_by(sp) %>% mutate(MedianValue = median(StdValue)) %>% ungroup()

hei.bubble <- master.new.superfinal %>% filter(TraitShort == "PlantHeight") %>%
  group_by(sp) %>% mutate(MedianValue = median(StdValue)) %>% ungroup()



## GAP-FILL PREP ----
sla.sp <- as.data.frame(unique(sla.bubble$sp)) #57 species
colnames(sla.sp) <- "SLA"
sla.sp.vector <- unique(sla.bubble$sp)

seed.sp <- as.data.frame(unique(seed.bubble$sp)) #28 species
colnames(seed.sp) <- "Seed_Mass"
seed.sp.vector <- unique(seed.bubble$sp)

hei.sp <- as.data.frame(unique(hei.bubble$sp)) #52 species
colnames(hei.sp) <- "Plant_Height"
hei.sp.vector <- unique(hei.bubble$sp)


## Which species have which trait data?
sla.seed2 <- seed.sp.vector %in% sla.sp.vector # "Thymus praecox" and "Thymus serpyllum" have seed mass but not SLA data.

hei.seed <- hei.sp[hei.sp$Plant_Height %in% seed.sp$Seed_Mass,]
hei.seed2 <- seed.sp.vector %in% hei.sp.vector # "Salix myrsinites" "Thymus serpyllum" and "Vaccinium oxycoccos" have seed mass but not height data. 



# We want to gap-fill those species with no seed data that have both height and SLA data
sla.hei <- hei.sp[hei.sp$Plant_Height %in% sla.sp$SLA,]
sla.hei2 <- hei.sp.vector %in% sla.sp.vector 
# From the 52 species with height data, all have SLA data except for:
# "Salix argyrocarpa", "Salix barrattiana" "Thymus praecox", "Vaccinium caespitosum" 

# 48 species that have both SLA and height data
sla.hei.df <- as.data.frame(unique(sla.hei))
colnames(sla.hei.df) <- "To_Fill"
write.csv(sla.hei.df, "data/2022_seed_data_needed.csv")



# identify the species that have SLA and height data but not seed mass 
hei.gap <- sla.hei.df %>% filter(To_Fill %notin% seed.sp.vector)
names(hei.gap)[names(hei.gap) == "To_Fill"] <- "sp"
hei.gap2 <- hei.gap %>% mutate(sp2 = sp) # 24 species

hei.gap.gen <- hei.gap2 %>% tidyr::separate(sp2, c("Genus", "SpeciesName"))
gap.genus.vector <- unique(hei.gap.gen$Genus) # 12 genera

# We want to gap-fill the list of species in "hei.gap" with seed data from "seed.sp"
# for each species in "hei.gap" we need to check if there are enough species in seed.sp to make the gap-fill possible



## SEED GAP-FILLING ----

# this is the file with trait data prior to merging
load("data/try5_clean5.RData") # [this file is not available in the repo as it contains raw data]

## creating this intermediate object just to inspect the gap-filling
gaps000 <- try5.clean5 %>% 
  tidyr::drop_na(Lat) %>% 
  filter(Lat > 60) %>% 
  filter(TraitShort == "SeedMass") %>%
  filter(Genus %in% gap.genus.vector) %>%
  group_by(Genus) %>%
  mutate(nobs = length(StdValue)) %>% 
  filter(nobs > 4)

# Check median and SD values per gap-filled genera
salixgaps <- gaps000 %>% filter(Genus == "Salix") %>% mutate(Median = median(StdValue)) %>% mutate(SD = sd(StdValue)) # 0.2188679, SD = 0.4349618
betugaps <- gaps000 %>% filter(Genus == "Betula") %>% mutate(Median = median(StdValue)) %>% mutate(SD = sd(StdValue)) #	0.11, SD = 0.2690477
vaccigaps <- gaps000 %>% filter(Genus == "Vaccinium") %>% mutate(Median = median(StdValue)) %>% mutate(SD = sd(StdValue)) # 0.22, SD = 0.1556895

# filtering for the genus and latitude we want, minimum of observations per genus 
# and calculate median and SD per genus
gaps <- try5.clean5 %>% 
  tidyr::drop_na(Lat) %>% 
  filter(Lat > 60) %>% 
  filter(TraitShort == "SeedMass") %>%
  mutate(LogTraitValue = log(StdValue)) %>%
  filter(Genus %in% gap.genus.vector) %>%
  group_by(Genus) %>%
  mutate(nobs = length(StdValue)) %>% 
  filter(nobs > 4) %>%
  mutate(GapMedian = median(LogTraitValue)) %>%
  mutate(GapSD = sd(LogTraitValue)) %>%
  dplyr::select(Genus, GapMedian, GapSD) %>%
  distinct(Genus, .keep_all = TRUE) %>%
  ungroup()

# Salix: median  -1.519287, SD 1.1615253
# Vaccinium: median -1.514128 SD 0.6117239
# Betula: median -2.211424, SD 0.9935412


# merge this with the species to gap-fill
seed.gap <- merge(hei.gap.gen, gaps, by = "Genus")

# create vector with species names
seed.vector <- unique(seed.gap$sp)

# save this file so we can merge it to the mastersheet
write.csv(seed.gap, "data/2022_gap_fill_seed_sp.csv")

# filter the rows we want so we leave one row only per species
master.seed <- master.new.superfinal %>% filter(sp %in% seed.vector) %>% 
  distinct(sp, .keep_all = TRUE) %>% mutate(Dataset = "Gap-filled") %>% 
  mutate(ValueKindName = "Gap-filled") %>% mutate(DataContributor = "Gap-filled") %>% 
  mutate(StdValue = NA) %>% mutate(Source = "Gap-filled") %>%
  mutate(nobs = NA) %>% ungroup() %>% mutate(TraitShort = "SeedMass")

# retain only the columns we need for merging so there are no issues  
gap.sp <- seed.gap %>% dplyr::select(sp, GapMedian, GapSD)

# merge the two datasets
gap.fill.seed <- merge(master.seed, gap.sp, by = "sp")

# save this file 
write.csv(gap.fill.seed, "data/2022_gapfilled_seed_data.csv")



## POLAR MAP ----

# combine cleaned trait databases into one
plant.data <- bind_rows(sla.bubble, hei.bubble, seed.bubble)

# make sure filter is correct
coord.data <- plant.data %>% filter(Lat > 30) %>% 
  tidyr::drop_na(Lon) %>% tidyr::drop_na(Lat)

# convert coords to utm
plant.coords <- transform_coord(x = coord.data, lon = "Lon", lat = "Lat",
                                new.names = c("lon.utm", "lat.utm"),
                                verbose = FALSE, bind = TRUE)

dt <- transform_coord(coord.data, bind = TRUE)


# plot map
(tundra.map <- basemap(data = dt, limits = 32, land.col = "#c9c7c1") + 
    geom_point(data = dt, aes(x = lon.proj, y = lat.proj, colour = TraitShort), 
               alpha = 0.4, size = 12) + labs(title = "(a)") +
    scale_colour_manual(values = c("#800080", "#E57E00","#0a6933"), 
                        guide = guide_legend(title = "Plant Traits", override.aes = list(alpha = 1)), 
                        labels = c("Height", "Seed Mass", "SLA")) +
    theme(legend.title = element_text(size = 24, face = "bold"), legend.text=element_text(size = 22),
          plot.title = element_text(face = "bold", size = 26)))






## BUBBLE PLOTS ----
gap.fill.seed <- read.csv("data/2022_gapfilled_seed_data.csv")

hei.sp2 <- hei.sp %>% mutate(SpeciesName = Plant_Height)
sla.sp2 <- sla.sp %>% mutate(SpeciesName = SLA)

all.traits.na <- full_join(hei.sp2, sla.sp2, by = "SpeciesName")

# Species with height data and no SLA data: Salix argyrocarpa, Salix barrattiana, Thymus praecox and Vaccinium caespitosum

# Species with SLA data and no Height data: Andromeda glaucophylla, Artemisia frigida, Artemisia gmelinii,
# Ribes glandulosum, Rosa acicularis, Rosa majalis, Salix myrsinites and Salix rosmarinifolia, Vaccinium oxycoccos

# Put together the two seed databases so we can compare them
seed.sp <- as.data.frame(unique(seed.bubble$sp)) #28 species
colnames(seed.sp) <- "Seed_Mass"

gap.sp.col <- as.data.frame(unique(gap.fill.seed$sp)) #12 species
colnames(gap.sp.col) <- "Seed_Mass"

all.seed.sp <- rbind(seed.sp, gap.sp.col)
all.seed.sp2 <- all.seed.sp %>% mutate(SpeciesName = Seed_Mass)

all.traits.na2 <- full_join(all.traits.na, all.seed.sp2, by = "SpeciesName")


# Species with no seed data but yes SLA data: Artemisia dracunculus 
# Chamaedaphne calyculata	Cornus canadensis	Cornus sericea	Dasiphora fruticosa	Diapensia lapponica Harrimanella hypnoides
# Rhododendron lapponicum Rhododendron tomentosum Ribes nigrum	Ribes spicatum Rubus chamaemorus
# Andromeda glaucophylla	Artemisia frigida Artemisia gmelinii Ribes glandulosum Rosa acicularis
# Rosa majalis Salix rosmarinifolia


# Species with no seed data but yes height data: Artemisia dracunculus, Chamaedaphne calyculata	
# Cornus canadensis	Cornus sericea	Dasiphora fruticosa	Diapensia lapponica Harrimanella hypnoides	
#  Rhododendron lapponicum Rhododendron tomentosum Ribes nigrum	Ribes spicatum	Rubus chamaemorus Salix argyrocarpa
#	Salix barrattiana	Vaccinium caespitosum


## SLA but no height
sla.nohei <- sla.bubble %>% 
  filter(sp %in% c("Andromeda glaucophylla", "Artemisia frigida",	"Artemisia gmelinii", 
                   "Ribes glandulosum", "Rosa acicularis", "Rosa majalis", "Salix rosmarinifolia", 
                   "Salix myrsinites", "Vaccinium oxycoccos")) %>% 
  distinct(sp, .keep_all = TRUE) %>% dplyr::select(sp, gf, MedianValue, StdValue) %>% 
  mutate(MedianValue = NA) %>% mutate(StdValue = NA)
  
hei.bubble2 <- hei.bubble %>% dplyr::select(sp, gf, MedianValue, StdValue)

# Add Thymus serpyllum which doesn't have either SLA or height but it has seed mass
thymus.serpyllum <- data.frame(sp = "Thymus serpyllum", gf = "Dwarf shrub", 
                               MedianValue = NA, StdValue = NA)

final.hei.bubble <- rbind(hei.bubble2, sla.nohei, thymus.serpyllum)
unique(final.hei.bubble$sp) #62 species


## Height but no SLA
hei.nosla <- hei.bubble %>% filter(sp %in% c("Salix argyrocarpa", "Salix barrattiana", 
                                             "Vaccinium caespitosum", "Thymus praecox")) %>% 
  distinct(sp, .keep_all = TRUE) %>% dplyr::select(sp, gf, MedianValue, StdValue) %>%  
  mutate(MedianValue = NA) %>% mutate(StdValue = NA)

sla.bubble2 <- sla.bubble %>% dplyr::select(sp, gf, MedianValue, StdValue) 
final.sla.bubble <- rbind(sla.bubble2, hei.nosla, thymus.serpyllum)
sort(unique(final.sla.bubble$sp)) # 62 species


# Put together the two seed databases so we can plot them
seed1 <- seed.bubble %>% dplyr::select(sp, gf, MedianValue, Source, StdValue)
seed2 <- gap.fill.seed %>% dplyr::select(sp, gf, GapMedian, Source) %>% 
  rename(MedianValue = GapMedian) %>% mutate(StdValue = NA) %>% 
  mutate(MedianValue = case_when(sp %in% c("Betula glandulosa", "Betula humilis") ~ 0.11,
                                 sp == "Vaccinium myrtilloides" ~ 0.22,
                                 TRUE ~ 0.2188679))
all.seeds <- rbind(seed1, seed2) # input the median (no-log) manually above just for bubble plots (values from vaccigaps, betugap, etc.)



## No seed but yes SLA and/or height
noseed <- master.new.superfinal %>% 
  filter(sp %in% c("Artemisia dracunculus", "Chamaedaphne calyculata",
                   "Cornus canadensis",	"Cornus sericea",	"Dasiphora fruticosa", 
                   "Diapensia lapponica", "Harrimanella hypnoides",
                   "Rhododendron tomentosum",	"Rhododendron lapponicum",
                   "Ribes nigrum",	"Ribes spicatum", "Rubus chamaemorus", "Andromeda glaucophylla",
                   "Artemisia frigida",	"Artemisia gmelinii", "Ribes glandulosum", "Rosa acicularis", "Rosa majalis",
                   "Salix rosmarinifolia", "Salix argyrocarpa", "Salix barrattiana", "Vaccinium caespitosum")) %>%
           distinct(sp, .keep_all = TRUE) %>% dplyr::select(sp, gf) %>%  
           mutate(MedianValue = NA) %>% mutate(Source = "No") %>% mutate(StdValue = NA)

final.seed.bubble <- rbind(all.seeds, noseed)
unique(final.seed.bubble$sp) # 62 species


# plotting all SLA values per species
(sla.bubbles.plot <- final.sla.bubble %>% 
    ggplot() +
    geom_point(aes(x = interaction(sp, gf), y = StdValue), size = 4, colour = "#0a6933", alpha = 0.7) + 
    geom_point(aes(x = interaction(sp, gf), y = MedianValue), size = 6, colour = "black") + 
    xlab("\nSpecies") + ylab(expression(bold(paste("SLA (mm"^bold("2"), "/mg)")))) + 
    scale_x_discrete("Species", breaks = interaction(final.sla.bubble$sp, final.sla.bubble$gf), 
                     labels = final.sla.bubble$sp) + ylim(0, 70) + 
    geom_vline(xintercept = 18.5, linetype = "solid") +
    geom_vline(xintercept = 28.5, linetype = "solid") +
    annotate("text", x = 10, y = 65, label = "Dwarf shrubs", size = 8) +
    annotate("text", x = 24, y = 65, label = "Low shrubs", size = 8) +
    annotate("text", x = 45, y = 65, label = "Tall shrubs", size = 8) +
    range.theme + 
    theme(axis.text.x  = element_text(face = "italic", angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text=element_text(size = 16),
          legend.position = c(0.85, 0.8), legend.key = element_blank(),
          legend.background = element_blank()))



# plotting all height values per species
(hei.bubbles.plot <- final.hei.bubble %>% 
   ggplot() +
   geom_point(aes(x = interaction(sp, gf), y = StdValue), size = 4, colour = "#800080", alpha = 0.5) + 
   geom_point(aes(x = interaction(sp, gf), y = MedianValue), size = 6, colour = "black") + 
   xlab("\nSpecies") + ylab("Plant Height (m)\n") + 
   scale_x_discrete("Species", breaks = interaction(final.hei.bubble$sp, final.hei.bubble$gf), labels = final.hei.bubble$sp) +
   geom_vline(xintercept = 18.5, linetype = "solid") +
   geom_vline(xintercept = 28.5, linetype = "solid") +
   annotate("text", x = 10, y = 3.8, label = "Dwarf shrubs", size = 8) +
   annotate("text", x = 24, y = 3.8, label = "Low shrubs", size = 8) +
   annotate("text", x = 45, y = 3.8, label = "Tall shrubs", size = 8) +
   range.theme + 
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text=element_text(size = 16),
         legend.position = c(0.85, 0.8), legend.key = element_blank(),
         legend.background = element_blank()))




# plotting all seed values per species
unique(final.seed.bubble$Source)
proper.sources <- c("LedaKleyer", "TTT", "TRY", "BobHollisterSeeds", "LucieSmrzova", "EstherLevesqueBylot", "No")

gap.points.bubble <- filter(final.seed.bubble, Source == "Gap-filled")
other.points.bubble <- filter(final.seed.bubble, Source %in% proper.sources)


(seed.bubbles.plot <-  
    ggplot() +
    geom_point(data = final.seed.bubble, aes(x = interaction(sp, gf), y = StdValue), size = 4, colour = "#E57E00", alpha = 0.7) + 
    geom_point(data = final.seed.bubble, aes(x = interaction(sp, gf), y = MedianValue), size = 6, colour = "black") + 
    geom_point(data = gap.points.bubble, aes(x = interaction(sp, gf), y = MedianValue), 
               pch = 21, fill="white", color="black", size = 5, stroke = 1.5) + 
    xlab("\nSpecies") + ylab("Seed Mass (mg)\n") + 
    scale_x_discrete("Species", breaks = interaction(final.seed.bubble$sp, final.seed.bubble$gf), labels = final.seed.bubble$sp) +
    geom_vline(xintercept = 18.5, linetype = "solid") +
    geom_vline(xintercept = 28.5, linetype = "solid") +
    annotate("text", x = 10, y = 23, label = "Dwarf shrubs", size = 8) +
    annotate("text", x = 24, y = 23, label = "Low shrubs", size = 8) +
    annotate("text", x = 45, y = 23, label = "Tall shrubs", size = 8) +
    range.theme + 
    theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
          legend.title = element_text(size = 18, face = "bold"), 
          legend.text=element_text(size = 16),
          legend.position = c(0.85, 0.8), legend.key = element_blank(),
          legend.background = element_blank()))




## PANEL (FIGURE 2) ----
(bubble.panel.cow <- plot_grid(hei.bubbles.plot, seed.bubbles.plot, sla.bubbles.plot, 
                            ncol=1, nrow = 3, align="v", axis = "bl", label_size = 26, 
                            labels = c("(b)", "(c)", "(d)")))

# I didn't manage to get the right proportions in R so I am putting together the map and the bubble plots in Photoshop.
