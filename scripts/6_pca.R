## Trait-range manuscript
## Mariana Garcia
## October 2019
## Script 6. Principal Component Analysis


#### LIBRARIES ----
library(tidyverse)
library(brms)
library(ggpubr)
library(vegan)
library(AMR) 
library(cowplot)
library(RVAideMemoire)


#### PREP DATA ----

# load species trait variation
trait.var.short <- read.csv("data/2022_all_sp_trait_var.csv")

# load species trait values
trait.val.short <- read.csv("data/2022_all_sp_trait_val.csv")

# select relevant columns only
trait.var <- trait.var.short %>% 
  dplyr::select(sp, Hei_var, Seed_var, SLA_var) %>% 
  rename(Height = Hei_var) %>% rename(SeedMass = Seed_var) %>% rename(SLA = SLA_var)

# select relevant columns only
trait.val <- trait.val.short %>% 
  dplyr::select(sp, Hei_val, SeedMass_val, SLA_val) %>%
  rename(Height = Hei_val) %>% rename(SeedMass = SeedMass_val) %>% rename(SLA = SLA_val)

# in alphabetical order - variation
trait.var.abc <- trait.var[order(trait.var$sp),]

# in alphabetical order - values
trait.val.abc <- trait.val[order(trait.val$sp),]

# convert species column into row names for plotting later
goodone <- trait.var.abc[-1]
row.names(goodone) <- trait.var.abc$sp
as.data.frame(goodone)

# Save species name in a vector
sp.names <- rownames(goodone)

# filter the functional groups & categories
fg <- trait.var.short %>% dplyr::select(sp, gf, category.abs)

# in alphabetical order
fg.abc <- fg[order(fg$sp),]

  


#### VARIATION PCA ----
traits.var.pca <- prcomp(trait.var.abc[,c(2:4)], center = TRUE, scale. = TRUE)

# scores of the observations ($x)
head(traits.var.pca)

# PC1 explains 49% of variance and PC2 explains 31% variance
summary(traits.var.pca)

# center point, scaling & SD
str(traits.var.pca)

# plot method returns a plot of the variances (y) associated with the PCs (x)
# useful to decide how many PCs to retain for further analysis
plot(traits.var.pca, type = "l")


## Plotting PCA
(pca.var.plot <- ggplot_pca(traits.var.pca, ellipse = TRUE, ellipse.prob = 0.68, scale = 1, labels = NULL, 
                                 groups = as.character(fg.abc$category), colour = as.character(fg.abc$category), 
                                 points_size = 2, arrows_colour = "black", labels_textsize = 0.001,
                                 arrows_size = 1, arrows_textsize = 10) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = fg.abc$category), type = "norm", level = 0.68) +
    geom_point(aes(colour = fg.abc$category), size = 10) +
    scale_colour_manual(name="Category", values = c("#CC6677", "#9eabad", "#44AA99")) +
    scale_fill_manual(name="Category", values = c("#CC6677", "#9eabad", "#44AA99")) +
    guides(color = guide_legend(override.aes = list(fill=NA))) +
    theme(axis.title.x = element_text(face="bold", size=28),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none"))




#### VALUES PCA ----
values.pca <- prcomp(trait.val.abc[,c(2:4)], center = TRUE, scale. = TRUE)
head(values.pca)

# PC1 explains 44% of variance and PC2 explains 30% variance
summary(values.pca)

# Plot PCA
(values.pca.plot.points <- ggplot_pca(values.pca, ellipse = TRUE, ellipse.prob = 0.68, scale = 1, labels = NULL, 
               groups = as.character(fg.abc$category), colour = as.character(fg.abc$category), 
               points_size = 2, arrows_colour = "black", labels_textsize = 0.001,
               arrows_size = 1, arrows_textsize = 10) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = fg.abc$category), type = "norm", level = 0.65) +
    geom_point(aes(colour = fg.abc$category), size = 10) +
    scale_colour_manual(name="Category", values = c("#CC6677", "#9eabad", "#44AA99")) +
    scale_fill_manual(name="Category", values = c("#CC6677", "#9eabad", "#44AA99")) +
    guides(color = guide_legend(override.aes = list(fill=NA))) +
    theme(axis.title.x = element_text(face="bold", size=28),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "bottom", legend.background = element_blank(), 
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 28), 
          legend.key = element_rect(fill = "transparent", colour = "transparent")))





#### PERMANOVA VARIATION ----
var_c <- scale(trait.var.abc[,c(2:4)], center = TRUE, scale = TRUE)

var.ad <- adonis2(var_c ~ category.abs, data = fg.abc, method = 'euclidean', permutations=999) # no significant difference here
print(as.data.frame(var.ad$aov.tab)["category.abs", "Pr(>F)"])

## Euclidean distances between samples (because we have negative values)
dis <- vegdist(var_c, "euclidean")
summary(dis)

## Calculate multivariate dispersions (variances; average distance to centroids)
mod <- betadisper(dis, fg.abc$category.abs)
mod
summary(mod)

# we have to use group dispersions to perform an ANOVA test. 
anova(mod) # p < 0.05 so differences between groups

# A Tukey's test can be done to see if and which groups differ in relation to their variances 
TukeyHSD(mod)
plot(mod) 

# Winner-No change groups are different
permutest(mod, pairwise = TRUE, permutations = 99)
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)




#### PERMANOVA VALUES ----
val_c <- scale(trait.val.abc[,c(2:4)], center = TRUE, scale = TRUE)

val.ad <- adonis2(val_c ~ category.abs, data = fg.abc, method = 'euclidean', permutations=999) # no significant difference
print(as.data.frame(val.ad$aov.tab)["category.abs", "Pr(>F)"])


## Euclidean distances between samples (because we have negative values)
dis.val <- vegdist(val_c, "euclidean")
summary(dis.val)

## Calculate multivariate dispersions (variances; average distance to centroids):
mod.val <- betadisper(dis.val, fg.abc$category.abs)
mod.val
summary(mod.val)

# We use group dispersions to perform an ANOVA test
anova(mod.val) # p > 0.05 so means homogeneity of variances between groups

# A Tukey's test can be done to see if and which groups differ in relation to 
# their variances (in our case, groups are not different):
TukeyHSD(mod.val)
plot(mod.val)


### Pairwise comparisons for groups for trait variation
var.win <- pairwise.perm.manova(var_c, fg.abc$category.abs, nperm=1000)
# no significant differences between the groups (p-value is > 0.05)

var.win2 <- pairwise.perm.manova(var_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                                 nperm = 1000, progress = TRUE, p.method = "bonferroni")
# no significant differences between the groups

var.win3 <- pairwise.perm.manova(var_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                                 nperm = 1000, progress = TRUE, p.method = "holm")
# no significant differences between the group




### Pairwise comparisons for groups for trait values
val.win <- pairwise.perm.manova(val_c, fg.abc$category.abs, nperm=999)
# ns

# Different p-correction methods
val.win3 <- pairwise.perm.manova(val_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                       nperm = 999, progress = TRUE, p.method = "fdr")
# ns

val.win4 <- pairwise.perm.manova(val_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                                 nperm = 999, progress = TRUE, p.method = "none")
# Marginally significant between winner and loser (p = 0.06)


val.win5 <- pairwise.perm.manova(val_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                                 nperm = 999, progress = TRUE, p.method = "bonferroni")
# ns

val.win6 <- pairwise.perm.manova(val_c, fg.abc$category.abs, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"),
                                 nperm = 999, progress = TRUE, p.method = "holm")
# ns



#### Figure 6 (panel)
(pca.panel.cow.full <- plot_grid(values.pca.plot.points, coefs.val.plot, pca.var.plot, 
                                 coefs.var.plot,
                            ncol=2, nrow=2, align="h", label_size = 26, 
                            labels = c("(a) PCA Trait values", "(c) Binomial model - Trait values", 
                                       "(b) PCA Trait variation", "(d) Binomial model - Trait variation")))

ggsave(pca.panel.cow.full, filename = "figures/Figure_6.png", width = 60, height = 60, units = "cm")






## VALUES MODELS ----

# PCA scores as model predictors (trait values) 

# Extract values per individual species
coords <- as.data.frame(values.pca$x)
fg.abc.t <- as.data.frame(fg.abc)

# Merge with original list of species
coord.table <- cbind(fg.abc.t, coords)

# List the species within the PCA in a vector
pca.vector.sp <- unique(coord.table$sp)

# load database with range values
load("data/2022_master_new_superfinal.RData")

# filter range data for our species
pca.sp <- master.new.superfinal %>% ungroup() %>%
  filter(sp %in% pca.vector.sp) %>% 
  dplyr::select(sp, CentredRangeLog, cent.rel.median.log, cent.abs.median.log) %>%
  distinct(sp, .keep_all = TRUE) 

# combine with PCA trait coords
lm.data <- merge(coord.table, pca.sp, by = "sp")



## Current range models
pca.lm <- brm(CentredRangeLog ~ PC1 + PC2, data = lm.data, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_current_pca")
summary(pca.lm) #  negative and positive

pca.lm2 <- brm(CentredRangeLog ~ PC1, data = lm.data, iter = 2000, chains = 4, warmup = 400, 
               file = "models/2022_current_pca1")
summary(pca.lm2) # negative ns




## Absolute range change models
abs.lm <- brm(cent.abs.median.log ~ PC1 + PC2, data = lm.data, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_abs_rc_pca")
summary(abs.lm) # both negative ns

abs.lm2 <- brm(cent.abs.median.log ~ PC1, data = lm.data, iter = 2000, chains = 4, warmup = 400, 
               file = "models/2022_abs_rc_pca2")
summary(abs.lm2) # negative ns




## Relative range change models
rel.lm <- brm(cent.rel.median.log ~ PC1 + PC2, data = lm.data, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_rel_rc_pca")
summary(rel.lm) # positive and negative ns

rel.lm2 <- brm(cent.rel.median.log ~ PC1, data = lm.data, iter = 2000, chains = 4, warmup = 400,
               file = "models/2022_rel_rc_pca2")
summary(rel.lm2) # positive ns




## Categorical models
pca.cat.mod <- brm(PC1 ~ category.abs, data = lm.data, iter = 2000, chains = 4, warmup = 400,
                   file = "models/2022_cat_pca1")
summary(pca.cat.mod)
conditional_effects(pca.cat.mod) #ns, overlap between all categories

pca.cat.mod2 <- brm(PC2 ~ category.abs, data = lm.data, iter = 2000, chains = 4, warmup = 400,
                    file = "models/2022_cat_pca2")
summary(pca.cat.mod2)
conditional_effects(pca.cat.mod2) #ns, overlap between all categories





## VARIATION MODELS ----

# Extract PCA values per individual species
coords.var <- as.data.frame(traits.var.pca$x)

# Merge with original list of species
coord.var.table <- cbind(fg.abc.t, coords.var)

# combine with range data
lm.var.data <- merge(coord.var.table, pca.sp, by = "sp")


## Current range models
pca.var.lm <- brm(CentredRangeLog ~ PC1 + PC2, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_current_pca_var")
summary(pca.var.lm) # negative ns


pca.var.lm2 <- brm(CentredRangeLog ~ PC1, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                  file = "models/2022_current_pca_var2")
summary(pca.var.lm2) # negative ns





## Absolute range change models
abs.var.lm <- brm(cent.abs.median.log ~ PC1 + PC2, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_abs_rc_pca_var")
summary(abs.var.lm) # positive and negative ns

abs.var.lm2 <- brm(cent.abs.median.log ~ PC1, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_abs_rc_pca_var2")
summary(abs.var.lm2) # positive ns



## Relative range change models
rel.var.lm <- brm(cent.rel.median.log ~ PC1 + PC2, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                  file = "models/2022_rel_rc_pca_var")
summary(rel.var.lm) # positive ns and negative significant

rel.var.lm2 <- brm(cent.rel.median.log ~ PC1, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                   file = "models/2022_rel_rc_pca_var2")
summary(rel.var.lm2) #ns




## Categorical models
pca.cat.var.mod <- brm(PC1 ~ category.abs, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                       file = "models/2022_cat_pca1_var")
summary(pca.cat.var.mod) #ns, overlap between all categories
conditional_effects(pca.cat.var.mod)

pca.cat.var.mod2 <- brm(PC2 ~ category.abs, data = lm.var.data, iter = 2000, chains = 4, warmup = 400, 
                        file = "models/2022_cat_pca2_var")
summary(pca.cat.var.mod2) #ns, overlap between all categories
conditional_effects(pca.cat.var.mod2)




#### COVER VALUES PCA ---- 
cov.slop <- read.csv("data/CoverChangeSpecies_QuantsCats_Final.csv")
names(cov.slop)[names(cov.slop) == 'Name'] <- 'sp'

# load species trait values
trait.val.short <- read.csv("data/2022_all_sp_trait_val.csv")

# merge trait values with cover category
cover.traits <- left_join(trait.val.short, cov.slop, by = "sp")
cover.traits2 <- tidyr::drop_na(cover.traits)

# select relevant columns only
cover.traits3 <- cover.traits2 %>% 
  dplyr::select(sp, Hei_val, SeedMass_val, SLA_val) %>%
  rename(Height = Hei_val) %>% rename(SeedMass = SeedMass_val) %>% rename(SLA = SLA_val)

# in alphabetical order - values
cover.traits3.abc <- cover.traits3[order(cover.traits3$sp),]

# convert species column into row names for plotting later
cover.name <- cover.traits3.abc[-1]
row.names(cover.name) <- cover.traits3.abc$sp
as.data.frame(cover.name)

# Save species name in a vector
sp.names <- rownames(cover.name)

# filter the functional groups & categories
fg <- cover.traits2 %>% dplyr::select(sp, gf, CoverCategory)

# in alphabetical order
fg.abc <- fg[order(fg$sp),]

# compute PCA values
cover.pca <- prcomp(cover.traits3.abc[,c(2:4)], center = TRUE, scale. = TRUE)
head(cover.pca)

# PC1 explains 41% of variance and PC2 explains 33 variance
summary(cover.pca)


## Plot points
(cover.pca.plot.points <- ggplot_pca(cover.pca, ellipse = TRUE, ellipse.prob = 0.68, scale = 1, labels = NULL, 
                                     groups = as.character(fg.abc$CoverCategory), colour = as.character(fg.abc$CoverCategory), 
                                     points_size = 2, arrows_colour = "black", labels_textsize = 0.001,
                                     arrows_size = 1, arrows_textsize = 9) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(colour = fg.abc$CoverCategory), size = 10) +
    stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = fg.abc$CoverCategory), type = "norm", level = 0.65) +
    scale_colour_manual(name="Category (cover change)", values = c("#CC6677", "#9eabad", "#44AA99")) +
    scale_fill_manual(name="Category (cover change)", values = c("#CC6677", "#9eabad", "#44AA99")) +
    guides(color = guide_legend(override.aes = list(fill=NA))) +
    theme(axis.title.x = element_text(face="bold", size=28),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "bottom", legend.background = element_blank(), 
          legend.title = element_text(size = 28, face = "bold"),
          legend.text = element_text(size = 28), 
          legend.key = element_rect(fill = "transparent", colour = "transparent")))





#### COVER VARIATION PCA -----

# load species trait variation
trait.var.short <- read.csv("data/2022_all_sp_trait_var.csv")

# merge trait values with cover category
cover.traits.var <- left_join(trait.var.short, cov.slop, by = "sp")
cover.traits.var2 <- tidyr::drop_na(cover.traits.var)

# select relevant columns only
cover.trait.vars3 <- cover.traits.var2 %>% 
  dplyr::select(sp, Hei_var, Seed_var, SLA_var) %>% 
  dplyr::rename(Height = Hei_var) %>% dplyr::rename(SeedMass = Seed_var) %>% dplyr::rename(SLA = SLA_var)


# in alphabetical order - variation
cover.trait.vars.abc <- cover.trait.vars3[order(cover.trait.vars3$sp),]

# convert species column into row names for plotting later
cover.var.name <- cover.trait.vars.abc [-1]
row.names(cover.var.name) <- cover.trait.vars.abc$sp
as.data.frame(cover.var.name)

# compute PCA variation
cover.var.pca <- prcomp(cover.trait.vars.abc[,c(2:4)], center = TRUE, scale. = TRUE)
head(cover.var.pca)

# PC1 explains 55% of variance and PC2 explains 29% variance
summary(cover.var.pca)

## Plot points
(cover.pca.var.plot <- ggplot_pca(cover.var.pca, ellipse = TRUE, ellipse.prob = 0.68, scale = 1, labels = NULL, 
                        groups = as.character(fg.abc$CoverCategory), colour = as.character(fg.abc$CoverCategory), 
                        points_size = 2, arrows_colour = "black", labels_textsize = 0.001,
                        arrows_size = 1, arrows_textsize = 9) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(colour = fg.abc$CoverCategory), size = 10) +
    stat_ellipse(geom = "polygon", alpha = 0.15, aes(fill = fg.abc$CoverCategory), type = "norm", level = 0.65) +
    scale_colour_manual(name="Category (cover change)", values = c("#CC6677", "#9eabad", "#44AA99")) +
    scale_fill_manual(name="Category (cover change)", values = c("#CC6677", "#9eabad", "#44AA99")) +
    guides(color = guide_legend(override.aes = list(fill=NA))) +
    theme(axis.title.x = element_text(face="bold", size=28),
          axis.text.x  = element_text(vjust=0.5, size=28, colour = "black"), 
          axis.title.y = element_text(face="bold", size=28),
          axis.text.y  = element_text(vjust=0.5, size=28, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), 
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "bottom", legend.background = element_blank(), 
          legend.title = element_text(size = 28, face = "bold"),
          legend.text = element_text(size = 28), 
          legend.key = element_rect(fill = "transparent", colour = "transparent")))


#### PANEL

# put PCA figures together
(pca.cover.panel <- plot_grid(cover.pca.plot.points + theme(legend.position="none"),
                              cover.pca.var.plot + theme(legend.position="none"),
                              ncol=2, align="h", label_size = 26, 
                              labels = c("(a) Trait values", "(b) Trait variation")))

# extract the legend from one of the plots
legendz <- get_legend(cover.pca.plot.points + 
                        guides(color = guide_legend(nrow = 1)) + 
                        theme(legend.position = "bottom") + 
                        guides(color=guide_legend(override.aes=list(fill=NA))))


# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
(pca.cover.panel2 <- plot_grid(pca.cover.panel, legendz, ncol = 1, rel_heights = c(1, .1)))

ggsave(pca.cover.panel2, filename = "figures/Figure_S5.png", 
       width = 60, height = 35, units = "cm")



#### PERMANOVA VARIATION COVER ----
var_c_cov <- scale(cover.trait.vars.abc[,c(2:4)], center = TRUE, scale = TRUE)

var.ad.cov <- adonis2(var_c_cov ~ CoverCategory, data = fg.abc, method = 'euclidean', permutations=999) # no significant difference here
print(as.data.frame(var.ad.cov$aov.tab)["CoverCategory", "Pr(>F)"])

## Euclidean distances between samples (because we have negative values)
dis_cov <- vegdist(var_c_cov, "euclidean")
summary(dis_cov)

## Calculate multivariate dispersions (variances; average distance to centroids):
mod_cov <- betadisper(dis_cov, fg.abc$CoverCategory)
mod_cov
summary(mod_cov)

# As mentioned before, we have to use group dispersions to perform an ANOVA test. 
anova(mod_cov) # p > 0.05 so means homogeneity of variances between groups

# A Tukey's test can be done to see if and which groups differ in relation to 
# their variances (in our case, groups are not different):
TukeyHSD(mod_cov)
plot(mod_cov)

permutest(mod_cov, pairwise = TRUE, permutations = 99)
(mod.HSD.cov <- TukeyHSD(mod_cov))
plot(mod.HSD.cov)



#### PERMANOVA VALUES COVER ----
val_c_cov <- scale(cover.traits3.abc[,c(2:4)], center = TRUE, scale = TRUE)

val.ad.cov <- adonis2(val_c_cov ~ CoverCategory, data = fg.abc, method = 'euclidean', permutations=999) # no sign diff
print(as.data.frame(val.ad.cov$aov.tab)["CoverCategory", "Pr(>F)"])


## Euclidean distances between samples (because we have negative values)
dis.val.cov <- vegdist(val_c_cov, "euclidean")
summary(dis.val.cov)

## Calculate multivariate dispersions (variances; average distance to centroids):
mod.val.cov <- betadisper(dis.val.cov, fg.abc$CoverCategory)
mod.val.cov
summary(mod.val.cov)

# As mentioned before, we have to use group dispersions to perform an ANOVA test. 
anova(mod.val.cov) # p > 0.05 so means homogeneity of variances between groups

# A Tukey's test can be done to see if and which groups differ in relation to 
# their variances (in our case, groups are not different):
TukeyHSD(mod.val.cov)
plot(mod.val.cov)



#### COVER VALUES MODELS ----
# PCA scores as model predictors (trait values) 

# Extract values per individual species
coords.cov <- as.data.frame(cover.pca$x)
fg.abc.t <- as.data.frame(fg.abc)

# Merge with original list of species
coord.table.cov <- cbind(fg.abc.t, coords.cov)

# List the species within the PCA in a vector
pca.vector.cov.sp <- unique(coord.table.cov$sp)

# combine with PCA trait coords
lm.data.cov <- merge(coord.table.cov, cov.slop, by = "sp")


## Cover change models
pca.cov.lm <- brm(MeanSlope ~ PC1 + PC2, data = lm.data.cov, iter = 2000, chains = 4, warmup = 400, 
              file = "models/2022_cov_pca_mod")
summary(pca.cov.lm) # ns


## Categorical models
pca.cat.cov.mod <- brm(PC1 ~ CoverCategory.x, data = lm.data.cov, iter = 2000, chains = 4, warmup = 400,
                   file = "models/2022_cov_cat_pca1")
summary(pca.cat.cov.mod) #ns
conditional_effects(pca.cat.cov.mod)

pca.cat.cov.mod2 <- brm(PC2 ~ CoverCategory.x, data = lm.data.cov, iter = 2000, chains = 4, warmup = 400,
                    file = "models/2022_cov_cat_pca2")
summary(pca.cat.cov.mod2) #ns
conditional_effects(pca.cat.cov.mod2)



#### COVER VARIATION MODELS ----

# Extract values per individual species
coords.var.cov <- as.data.frame(cover.var.pca$x)

# Merge with original list of species
coord.var.table.cov <- cbind(fg.abc.t, coords.var.cov)

# combine with PCA trait coords
lm.data.var.cov <- merge(coord.var.table.cov, cov.slop, by = "sp")


## Cover change models
pca.cov.var.lm <- brm(MeanSlope ~ PC1 + PC2, data = lm.data.var.cov, iter = 2000, chains = 4, warmup = 400, 
                  file = "models/2022_cov_pca_var_mod")
summary(pca.cov.lm) # ns


## Categorical models
pca.cat.cov.var.mod <- brm(PC1 ~ CoverCategory.x, data = lm.data.var.cov, iter = 2000, chains = 4, warmup = 400,
                       file = "models/2022_cov_cat_var_pca1")
summary(pca.cat.cov.var.mod) #ns, overlap between all categories
conditional_effects(pca.cat.cov.var.mod)

pca.cat.cov.var.mod2 <- brm(PC2 ~ CoverCategory.x, data = lm.data.var.cov, iter = 2000, chains = 4, warmup = 400,
                        file = "models/2022_cov_cat_var_pca2")
summary(pca.cat.cov.var.mod2) #ns, overlap between all categories
conditional_effects(pca.cat.cov.var.mod2)


