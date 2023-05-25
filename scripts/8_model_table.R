## Trait-range manuscript
## Mariana Garcia
## July 2020
## Script 8. Model table


## LIBRARIES ----
library(dplyr)
library(broom)
library(stargazer)


# Let's break this down since otherwise R crashes
# Load a few models in their respective scripts, run the table, and then close R and clean the environment
# Rinse and repeat with the rest of models

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}



## SCRIPTS 3&4 MODEL TABLE ----

## Add model objects to list
models.list34 <- list(
  sla.var.mod, berry.seed.var.mod, wind.seed.var.mod, seed.var.nogap.mod,
  seed.var.mod, height.var.mod,
  all.traits.var.mod, all.traits.var.mod3,
  sla.val.mod, seed.val.nogap.mod, berry.seed.mod, wind.seed.mod, 
  seed.val.mod, hei.val.mod,
  all.traits.val.mod, all.traits.val.mod4,
  hei.sla, hei.seed, seed.sla,
  sla.var.cat.mod, seed.var.cat.mod, hei.var.cat.mod,
  sla.val.cat.mod, seed.val.cat.mod, hei.val.cat.mod,
  current.fg.mod, current.disp.mod, current.decid.mod, current.fam.mod,
  all.traits.var.mod25, all.traits.val.mod234,
  all.traits.var.cov.mod, all.traits.val.cov.mod,
  ranges.brm, fut.cur.brm)
  
## Create dataframe with model numbers/names
model_names34 <- c(
  "Current ranges - SLA variation", "Current ranges - Seed variation Berry Sp", "Current ranges - Seed variation Wind Sp", "Current ranges - Seed variation No Gapfilled",
  "Current ranges - Seed Mass variation", "Current ranges - Height variation",
  "Current ranges - all traits variation", "Current ranges - all traits variation interaction", 
  "Current ranges - SLA values", "Current ranges - Seed values No Gapfilled", "Current ranges - Seed values Berry Sp", "Current ranges - Seed values Wind Sp",
  "Current ranges - Seed Mass values", "Current ranges - Height values",
  "Current ranges - all traits values", "Current ranges - all traits values interaction", 
  "Height/SLA correlation", "Height/Seed Mass correlation", "Seed Mass/SLA correlation",
  "SLA variation - category", "Seed Mass variation - category", "Height variation - category",
  "SLA values - category", "Seed Mass values - category", "Height values - category",
  "Current ranges - functional group", "Current ranges - dispersal", "Current ranges - decidiousness", "Current ranges - family",
  "Binomial - category - all traits variation", "Binomial - category - all traits values",
  "Binomial - cover category - all traits variation", "Binomial - cover category - all traits values",
  "Absolute vs relative range change", "Future vs current ranges")

# Bind them together
model_number34 <- 1:35
mod.df34 <- data.frame(model_number34, model_names34)

# Extract parameters
mod.table34 <- lapply(models.list34, p_summarize) %>% 
  bind_rows(.id = "model_number34") 

# Add model name to table
mod.table34$model_number34 <- as.integer(mod.table34$model_number34)
mod.table.final34 <- left_join(mod.table34, mod.df34, by = "model_number34")

# Clean model parameters
mod.table.final34.2 <- mod.table.final34 %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final34.2$model_names34[duplicated(mod.table.final34.2$model_names34)] <- "  "
mod.table.final34.2$model_names34 <- as.character(mod.table.final34.2$model_names34)
mod.table.final34.2$model_number34[duplicated(mod.table.final34.2$model_number34)] <- "  "

colnames(mod.table.final34.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final34.3 <- mod.table.final34.2[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final34.4 <- mod.table.final34.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final34.4, "models/table_s2/2022_table_s2_34.csv")

# Convert to table
stargazer(mod.table.final34.4, type = "html", summary = FALSE)





## SCRIPT 5, PART 1 ----
models.list51 <- list(
  sla.abs25.mod, sla.medabs.mod, sla.abs75.mod,
  hei.abs25.mod, hei.medabs.mod, hei.abs75.mod,
  seed.abs25.mod, berry.seed.var.rc.mod, wind.seed.var.rc.mod, seed.var.nogap.rc.mod,
  seed.medabs.mod, seed.abs75.mod,
  all.traits.var.abs.mod, all.traits.var.abs.mod4,
  
  sla.rel25.mod, sla.medrel.mod, sla.rel75.mod,
  hei.rel25.mod, hei.medrel.mod, hei.rel75.mod,
  seed.rel25.mod, seed.medrel.mod, seed.rel75.mod,
  all.traits.var.rel.mod, all.traits.var.rel.mod4,
  
  sla.val.abs25.mod, sla.val.medabs.mod, sla.val.75abs.mod,
  seed.val.25abs.mod, berry.seed.val.rc.mod, wind.seed.val.rc.mod, seed.val.nogap.rc.mod,
  seed.val.medabs.mod, seed.val.75abs.mod,
  hei.val.abs25.mod, hei.val.absmed.mod, hei.val.abs75.mod,
  all.traits.val.abs.mod, all.traits.val.abs.mod4,
  
  sla.val.rel25.mod, sla.val.relmed.mod, sla.val.rel75.mod,
  seed.val.rel25.mod, seed.val.relmed.mod, seed.val.rel75.mod)


model_names51 <- c(
  "SLA variation - 25% absolute range change", "SLA variation - median absolute range change", "SLA variation - 75% absolute range change",
  "Height variation - 25% absolute range change", "Height variation - median absolute range change", "Height variation - 75% absolute range change",
  "Seed Mass variation - 25% absolute range change", "Seed Mass variation - Median RCh Berry", "Seed Mass variation - Median RCh Wind", "Seed Mass variation - Median RCh NoGapfill",
  "Seed Mass variation - median absolute range change", "Seed Mass variation - 75% absolute range change",
  "Median absolute range change - all traits variation", "Median absolute range change - all traits variation interaction",
  
  "SLA variation - 25% relative range change", "SLA variation - median relative range change", "SLA variation - 75% relative range change",
  "Height variation - 25% relative range change", "Height variation - median relative range change", "Height variation - 75% relative range change",
  "Seed Mass variation - 25% relative range change", "Seed Mass variation - median relative range change", "Seed Mass variation - 75% relative range change",
  "Median relative range change - all traits variation", "Median relative range change - all traits variation interaction",
  
  "SLA values - 25% absolute range change", "SLA values - median absolute range change", "SLA values - 75% absolute range change",
  "Seed Mass values - 25% absolute range change", "Seed Mass values - Median RCh Berry", "Seed Mass values - Median RCh Wind", "Seed Mass values - Median RCh NoGapfill",
  "Seed Mass values - median absolute range change", "Seed Mass values - 75% absolute range change",
  "Height values - 25% absolute range change", "Height values - median absolute range change", "Height values - 75% absolute range change",
  "Median absolute range change - all traits values", "Median absolute range change - all traits values interaction",
  
  "SLA values - 25% relative range change", "SLA values - median relative  range change", "SLA values - 75% relative range change",
  "Seed Mass values - 25% relative range change", "Seed Mass values - median relative  range change", "Seed Mass values - 75% relative range change")


# Bind them together
model_number51 <- 36:80
mod.df51 <- data.frame(model_number51, model_names51)

# Extract model parameters
mod.table51 <- lapply(models.list51, p_summarize) %>% 
  bind_rows(.id = "model_number51")

# Add model name to table
mod.table51$model_number51 <- as.integer(mod.table51$model_number51)
mod.table51$model_number51 <- (mod.table51$model_number51 + 35) # to have the correct number
mod.table.final51 <- left_join(mod.table51, mod.df51, by = "model_number51")

# Clean model parameters
mod.table.final51.2 <- mod.table.final51 %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final51.2$model_names51[duplicated(mod.table.final51.2$model_names51)] <- "  "
mod.table.final51.2$model_names51 <- as.character(mod.table.final51.2$model_names51)
mod.table.final51.2$model_number51[duplicated(mod.table.final51.2$model_number51)] <- "  "

colnames(mod.table.final51.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final51.3 <- mod.table.final51.2[, c(1, 7, 2, 3, 4, 5, 6)]

# Round to 3 decimals only because not working on stargazer function
mod.table.final51.4 <- mod.table.final51.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final51.4, "models/table_s2/2022_table_s2_51.csv")

# Convert to table
stargazer(mod.table.final51.4, type = "html", summary = FALSE)








## SCRIPT 5 - PART 2 ----
models.list52 <- list(hei.val.rel25.mod, hei.val.relmed.mod, hei.val.rel75.mod,
                    all.traits.val.rel.mod, all.traits.val.rel.mod4,
                    
                    sla.var.gain.mod, sla.var.loss.mod,
                    seed.var.gain.mod, seed.var.loss.mod,
                    hei.var.gain.mod, hei.var.loss.mod,
                    all.traits.var.abs.gains.mod, all.traits.var.abs.gains.mod3,
                    all.traits.var.abs.losses.mod, all.traits.var.abs.losses.mod3,
                    
                    sla.val.gain.mod, sla.val.loss.mod,
                    seed.val.gain.mod, seed.val.loss.mod,
                    hei.val.gain.mod, hei.val.loss.mod, 
                    all.traits.val.abs.gains.mod, all.traits.val.abs.gains.mod3,
                    all.traits.val.abs.losses.mod, all.traits.val.abs.losses.mod3,
                    
                    sla.var.rel.gain.mod, sla.var.rel.loss.mod,
                    seed.var.rel.gain.mod, seed.var.rel.loss.mod,
                    hei.var.rel.gain.mod, hei.var.rel.loss.mod,
                    all.traits.var.rel.gains.mod, all.traits.var.rel.gains.mod3,
                    all.traits.var.rel.losses.mod, all.traits.var.rel.losses.mod3)


model_names52 <- c(
  "Height values - 25% relative range change", "Height values - median relative  range change", "Height values - 75% relative range change",
  "Median relative range change - all traits values", "Median relative range change - all traits values interaction",
  
  "Absolute range gains - SLA variation", "Absolute range losses - SLA variation",
  "Absolute range gains - Seed Mass variation", "Absolute range losses - Seed Mass variation",
  "Absolute range gains - Height variation", "Absolute range losses - Height variation",
  "Median absolute range gains - all traits variation", "Median absolute range gains - all traits variation interaction",
  "Median absolute range losses - all traits variation", "Median absolute range losses - all traits variation interaction",
  
  "Absolute range gains - SLA values", "Absolute range losses - SLA values",
  "Absolute range gains - Seed Mass values", "Absolute range losses - Seed Mass values",
  "Absolute range gains - Height values", "Absolute range losses - Height values",
  "Median absolute range gains - all traits values", "Median absolute range gains - all traits values interaction",
  "Median absolute range losses - all traits values", "Median absolute range losses - all traits values interaction",
  
  "Relative range gains - SLA variation", "Relative range losses - SLA variation",
  "Relative range gains - Seed Mass variation", "Relative range losses - Seed Mass variation",
  "Relative range gains - Height variation", "Relative range losses - Height variation",
  "Median relative range gains - all traits variation", "Median relative range gains - all traits variation interaction",
  "Median relative range losses - all traits variation", "Median relative range losses - all traits variation interaction")
  

# Bind them together
model_number52 <- 81:115
mod.df52 <- data.frame(model_number52, model_names52)

# Extract model parameters
mod.table52 <- lapply(models.list52, p_summarize) %>% 
  bind_rows(.id = "model_number52")

# Add model name to table
mod.table52$model_number52 <- as.integer(mod.table52$model_number52)
mod.table52$model_number52 <- (mod.table52$model_number52 + 80) # to have the correct number
mod.table.final52 <- left_join(mod.table52, mod.df52, by = "model_number52")

# Clean model parameters
mod.table.final52.2 <- mod.table.final52 %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final52.2$model_names52[duplicated(mod.table.final52.2$model_names52)] <- "  "
mod.table.final52.2$model_names52 <- as.character(mod.table.final52.2$model_names52)
mod.table.final52.2$model_number52[duplicated(mod.table.final52.2$model_number52)] <- "  "

colnames(mod.table.final52.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final52.3 <- mod.table.final52.2[, c(1, 7, 2, 3, 4, 5, 6)]

# Round to 3 decimals only because not working on stargazer function
mod.table.final52.4 <- mod.table.final52.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final52.4, "models/table_s2/2022_table_s2_52.csv")

# Convert to table
stargazer(mod.table.final52.4, type = "html", summary = FALSE)




## SCRIPTS 5, 6 & 7 ----

models.list567 <- list(
sla.val.rel.gain.mod, sla.val.rel.loss.mod,
seed.val.rel.gain.mod, seed.val.rel.loss.mod,
hei.val.rel.gain.mod, hei.val.rel.loss.mod,
all.traits.val.rel.gains.mod, all.traits.val.rel.gains.mod3,
all.traits.val.rel.loss.mod, all.traits.val.rel.loss.mod3,

future.fg.mod, future.disp.mod,
future.decid.mod, future.fam.mod,

rel.fg.mod, rel.disp.mod,
rel.decid.mod, rel.fam.mod,

sla.cov.var.mod, seed.cov.var.mod, hei.cov.var.mod,
all.traits.var.cov.mod, 
sla.cov.val.mod, seed.cov.val.mod, hei.cov.val.mod,
all.traits.val.cov.mod,

pca.lm, pca.lm2,
abs.lm, abs.lm2, 
rel.lm, rel.lm2, 
pca.cat.mod, pca.cat.mod2,

pca.var.lm, pca.var.lm2,
abs.var.lm, abs.var.lm2,
rel.var.lm, rel.var.lm2,
pca.cat.var.mod, pca.cat.var.mod2,

pca.cov.lm, pca.cat.cov.mod, pca.cat.cov.mod2,
pca.cov.var.lm, pca.cat.cov.var.mod, pca.cat.cov.var.mod2,

maxdisp.seed.sp.mod, unlim.seed.sp.mod,
maxdisp.hei.sp.mod, unlim.hei.sp.mod,

unlim.hei.sp.wei.mod, unlim.sla.sp.wei.mod, unlim.seed.sp.wei.mod,
hei.var.unlim.mod, sla.unlimvar.wei.mod, seed.unlimvar.mod
)



model_names567 <- c(
"Relative range gains - SLA values", "Relative range losses - SLA values",
"Relative range gains - Seed Mass values", "Relative range losses - Seed Mass values",
"Relative range gains - Height values", "Relative range losses - Height values",
"Median relative range gains - all traits values", "Median relative range gains - all traits values interaction",
"Median relative range losses - all traits values", "Median relative range losses - all traits values interaction",

"Absolute range changes - functional group", "Absolute range changes - dispersal", 
"Absolute range changes - decidiousness", "Absolute range changes - family",

"Relative range changes - functional group", "Relative range changes - dispersal", 
"Relative range changes - decidiousness", "Relative range changes - family",

"Cover change slopes - SLA variation", "Cover change slopes - Seed Mass variation", "Cover change slopes - Height variation",
"Cover change slopes - all traits variation",
"Cover change slopes - SLA values", "Cover change slopes - Seed Mass values", "Cover change slopes - Height values",
"Cover change slopes - all traits values",

"Current range - PC1+PC2 values", "Current range - PC1 values",
"Absolute range change - PC1+PC2 values", "Absolute range change - PC1 values",
"Relative range change - PC1+PC2 values", "Relative range change - PC1 values",
"PC1 values - Category", "PC2 values - Category",

"Current range - PC1+PC2 variation", "Current range - PC1 variation",
"Absolute range change - PC1+PC2 variation", "Absolute range change - PC1 variation",
"Relative range change - PC1+PC2 variation", "Relative range change - PC1 variation",
"PC1 variation - Category", "PC2 variation - Category",

"Cover change - PC1+PC2 values", "PC1 values - Cover Category", "PC2 values - Cover Category",
"Cover change - PC1+PC2 variation", "PC1 variation - Cover Category", "PC2 variation - Cover Category",

"Median absolute range shift - Seed variation LimDisp", "Median absolute range shift - Seed variation UnlimDisp",
"Median absolute range shift - Height variation LimDisp", "Median absolute range shift - Height variation UnlimDisp",

"Median absolute range shift - Height values UnlimDispWeight", "Median absolute range shift - SLA values UnlimDispWeight", "Median absolute range shift - SeedMass values UnlimDispWeight",
"Median absolute range shift - Height variation UnlimDispWeight", "Median absolute range shift - SLA variation UnlimDispWeight", "Median absolute range shift - SeedMass variation UnlimDispWeight")




# Bind them together
model_number567 <- 116:173
mod.df567 <- data.frame(model_number567, model_names567)

# Extract model parameters
mod.table567 <- lapply(models.list567, p_summarize) %>% 
  bind_rows(.id = "model_number567")

# Add model name to table
mod.table567$model_number567 <- as.integer(mod.table567$model_number567)
mod.table567$model_number567 <- (mod.table567$model_number567 + 115) # to have the correct number
mod.table.final567 <- left_join(mod.table567, mod.df567, by = "model_number567")

# Clean model parameters
mod.table.final567.2 <- mod.table.final567 %>% filter(parameter != "lp__") %>% filter(parameter != "Intercept")
mod.table.final567.2$model_names567[duplicated(mod.table.final567.2$model_names567)] <- "  "
mod.table.final567.2$model_names567 <- as.character(mod.table.final567.2$model_names567)
mod.table.final567.2$model_number567[duplicated(mod.table.final567.2$model_number567)] <- "  "

colnames(mod.table.final567.2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final567.3 <- mod.table.final567.2[, c(1, 7, 2, 3, 4, 5, 6)]

# Round to 3 decimals only because not working on stargazer function
mod.table.final567.4 <- mod.table.final567.3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final567.4, "models/table_s2/2022_table_s2_567.csv")

# Convert to table
stargazer(mod.table.final567.4, type = "html", summary = FALSE)

