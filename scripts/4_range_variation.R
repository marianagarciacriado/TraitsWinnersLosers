## Trait-range manuscript
## Mariana Garcia
## September 2019
## Script 4. Winners vs losers (range variability)


# The input files for this script are not available in the repo given that they include raw data.
# The summarised data is provided below as range_quan.csv in the repo and can be loaded to create the figures
# and models of this script.


## LIBRARIES ----
library(tidyverse)
library(ggpubr)
library(brms)
library(stargazer)



## THEME ----
range.theme2 <- theme(legend.position = "none",
                     axis.title.x = element_text(face="bold", size=22),
                     axis.text.x  = element_text(vjust=0.5, size=16, colour = "black"), 
                     axis.title.y = element_text(face="bold", size=22),
                     axis.text.y  = element_text(vjust=0.5, size=20, colour = "black"),
                     panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
                     panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                     plot.title = element_text(color = "black", size = 18, face = "bold", hjust = 0.5),
                     plot.margin = unit(c(1,1,1,1), units = , "cm"))



## DATA LOADING & PREP ----

# Import database [file not available as it contains raw data]
range.quan <- read.csv("data/2022_range_full.csv")

# Re-calculate the quantiles with the raw data 
range.quan2 <- range.quan %>%
    group_by(sp) %>% 
    mutate(rel.quan25 = quantile(SpeciesRangeChange, probs = 0.25)) %>%
    mutate(rel.median = quantile(SpeciesRangeChange, probs = 0.50)) %>%
    mutate(rel.quan75 = quantile(SpeciesRangeChange, probs = 0.75)) %>%
    mutate(abs.quan25 = quantile(AbsoluteRangeChangeKm, probs = 0.25)) %>%
    mutate(abs.median = quantile(AbsoluteRangeChangeKm, probs = 0.50)) %>%
    mutate(abs.quan75 = quantile(AbsoluteRangeChangeKm, probs = 0.75)) %>%
    mutate(future.med = quantile(FutureRangeSize.FullDisp.Km, probs = 0.50)) %>%
    distinct(sp, .keep_all = TRUE) %>% ungroup()


# Add the 'winner' and 'loser' categories for relative and absolute change
range.quan2$category.rel[range.quan2$rel.quan25 > 0] <- "Winner"
range.quan2$category.rel[range.quan2$rel.quan75 < 0] <- "Loser"
range.quan2$category.rel[range.quan2$rel.quan25 <= 0 & range.quan2$rel.quan75 >= 0] <- "No change"
range.quan2$category.abs[range.quan2$abs.quan25 > 0] <- "Winner"
range.quan2$category.abs[range.quan2$abs.quan75 < 0] <- "Loser"
range.quan2$category.abs[range.quan2$abs.quan25 <= 0 & range.quan2$abs.quan75 >= 0] <- "No change"


# Merge with correct functional groups
load("data/2022_trait_cat.RData")

# retain only species and FG
fg.data <- trait.cat %>% dplyr::select(sp, gf) %>% distinct(sp, .keep_all = TRUE)

# remove original FG column so there are no issues when merging
range.quan3 <- range.quan2 %>% dplyr::select(., -gf)

# merge dataframes
range.quan4 <- left_join(range.quan3, fg.data, by = "sp")

# Save this file 
write.csv(range.quan4, file = "data/range_quan.csv")




## FIGURES PER GROUP ----
range.quan4 <- read.csv("data/range_quan.csv")

# extract figures per category
stats.cat <- range.quan4 %>% group_by(category.abs) %>% summarise(n = n())

# extract figures per category and functional group
stats.cat2 <- range.quan4 %>% group_by(category.abs, gf) %>% summarise(n = n())

# let's see species per category
winns <- range.quan4 %>% filter(category.abs == "Winner")
losers <- range.quan4 %>% filter(category.abs == "Loser")
noch <- range.quan4 %>% filter(category.abs == "No change")



## TABLE FOR APPENDIX ----

# Select only the columns we need
wins.table <- range.quan4 %>% select(sp, abs.median, rel.median, category.abs, gf)  
wins.table2 <- wins.table[order(-wins.table$abs.median),]
str(wins.table2)
wins.table2$gf <- as.character(wins.table2$gf)

# Create nice table in html
stargazer(wins.table2, title = "Winners vs losers", type = "html", summary = FALSE, 
          out = "figures/Table_S1.htm")

# Make it usable in Word 
# Copy the html code that appears on the console, click on File/New/R HTML file
# Delete everything in the newly generated file except the <html> and </html> tags and include the copied html code there
# Click Knit, save the html file
# The file can now be opened with Word for final touches




## CATERPILLAR PLOT (RELATIVE RANGE CHANGE) ----

# first extract the 'no changers'
nochange.rel <- range.quan4 %>% filter(category.rel == "No change")


# plot this
(ranges.rel <- ggplot(range.quan4, aes(x = reorder(sp, -abs.median), y = rel.median, colour = gf)) + 
        geom_point(size = 6) + 
        geom_errorbar(aes(ymin = rel.quan25, ymax = rel.quan75), size = 1) + 
        geom_rect(data = nochange.rel, aes(xmin = 19.5, xmax = 20.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 22.5, xmax = 23.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 26.5, xmax = 27.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 29.5, xmax = 30.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 38.5, xmax = 39.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 40.5, xmax = 41.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 42.5, xmax = 43.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.rel, aes(xmin = 44.5, xmax = 45.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_hline(yintercept = 0, linetype = "solid") +
    annotate("text", x = 10, y = 230, label = "Winners", size = 8) +
    annotate("text", x = 55, y = 230, label = "Losers", size = 8) +
    xlab("\nSpecies") + ylab("Relative Species Range Change (%)\n") + 
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) + 
    scale_y_continuous(breaks = seq(-100, 240, 40)) + 
    labs(colour = "Functional group") + range.theme2 + 
    theme(axis.text.x  = element_text(face = "italic", angle = 45, 
                                      vjust = 1, hjust = 1,
                                      size = 16, colour = "black"), 
          legend.title = element_text(size = 25, face = "bold"), 
          legend.text=element_text(size = 25),
          legend.position = c(0.9, 0.7), legend.key = element_blank(),
          legend.background = element_blank()))




## CATERPILLAR PLOT (ABSOLUTE RANGE CHANGE) ----

# first extract the 'no changers' (they are the same in relative vs absolute)
nochange.abs <- range.quan4 %>% filter(category.abs == "No change")

# now plot the caterpillar graph
(ranges.abs <- ggplot(range.quan4, aes(x = reorder(sp, -abs.median), y = abs.median, colour = gf)) + 
    geom_point(size = 6) + 
    geom_errorbar(aes(ymin = abs.quan25, ymax = abs.quan75), size = 1) + 
        geom_rect(data = nochange.abs, aes(xmin = 19.5, xmax = 20.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 22.5, xmax = 23.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 26.5, xmax = 27.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 29.5, xmax = 30.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 38.5, xmax = 39.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 40.5, xmax = 41.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 42.5, xmax = 43.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
        geom_rect(data = nochange.abs, aes(xmin = 44.5, xmax = 45.5, ymin=-Inf, ymax=+Inf), 
                  col = 'white', fill = 'grey', alpha = 0.05) + 
    geom_hline(yintercept = 0, linetype = "solid") +
    annotate("text", x = 10, y = 17000000, label = "Winners", size = 8) +
    annotate("text", x = 55, y = 17000000, label = "Losers", size = 8) +
    xlab("\nSpecies") + ylab(expression(bold(paste("Absolute Species Range Change (km"^bold("2"), ")")))) +
    scale_colour_manual(values = c("#d2a5d2","#993299", "#4C004C")) + 
    scale_y_continuous(breaks = seq(-10000000, 18000000, 4000000),
                       labels = function(x) format(x, scientific = FALSE)) +
    labs(colour = "Functional group") + range.theme2 + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_text(face="bold", size=22),
          legend.title = element_text(size = 25, face = "bold"), 
          legend.text=element_text(size = 25),
          legend.position = c(0.9, 0.7), legend.key = element_blank(),
          legend.background = element_blank()))




## CATERPILLAR PANEL ----
(range.panel <- ggarrange(ranges.abs, ranges.rel,
                            labels = c("(a)", "(b)"), common.legend = TRUE, align = "v",
                          nrow = 2, ncol = 1, font.label = list(size = 26)))

ggsave(range.panel, filename = "figures/Figure_4.png", 
       width = 60, height = 60, units = "cm")




## ABSOLUTE RANGE CHANGE VS RELATIVE RANGE CHANGE ----

# Positive relationship between relative and absolute range changes
(abs.rel.plot <- ggplot(range.quan4, aes(x = rel.median, y = abs.median)) +
    geom_point(size = 4, colour = "blue") + 
    xlab("\nRelative range change (%)") + ylab("Absolute range change (sq km)\n") + 
    geom_smooth(method = 'lm', formula = y~x) + range.theme2)

## Fit model
ranges.brm <- brm(abs.median ~ rel.median, data = range.quan4, iter = 2000, chains = 4, warmup = 400, 
                  file = "models/rel_vs_abs_ranges")
summary(ranges.brm) # significant positive relationship


## FUTURE RANGE VS CURRENT RANGE ----

# [file range.quan2 is not available in the repo as it includes raw data]
(fut.cur.plot <- ggplot(range.quan2, aes(x = CurrentRangeSizeKm, y = future.med)) +
   geom_point(size = 4, colour = "turquoise") + 
   xlab("\nCurrent range size (sq km)") + ylab("Future range size (sq km)\n") + 
   geom_smooth(method = 'lm', formula = y~x, colour = "turquoise") + 
   scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
   scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) +
   range.theme2)

# Fit model
fut.cur.brm <- brm(future.med ~ CurrentRangeSize, data = range.quan2, iter = 2000, chains = 4, warmup = 400,
                   file = "models/fut_vs_cur_ranges")
summary(fut.cur.brm) # significant positive relationship
