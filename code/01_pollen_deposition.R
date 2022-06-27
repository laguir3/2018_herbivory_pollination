#### DESCRIPTION ####
# Analysis for Aguirre et al. Effects of herbivory on community-wide pollination processes.
# 
# Code by Luis A. Aguirre
#
#### LOAD ####
# Packages 
library(glmmTMB) 
library(tidyverse)
library(DHARMa)
library(emmeans)
library(car)
library(lmtest)

# Set color-blind palette
cb <- c("#000000", # black
        "#E69F00", # orange
        "#56B4E9", # light blue
        "#009E73", # green
        "#F0E442", # yellow
        "#0072B2", # dark blue
        "#D55E00", # red
        "#CC79A7") # pink

# Load datasets 
# Vetch
vetch18 <- read.csv("data/2018_vetch_deposition.csv", 
                    header = TRUE)

vetch19 <- read.csv("data/2019_vetch_deposition.csv", 
                    header = TRUE)

# Galium 
galium18 <- read.csv("data/2018_galium_deposition.csv", 
                     header = TRUE)

galium19 <- read.csv("data/2019_galium_deposition.csv", 
                     header = TRUE)


# Wild Basil
basil18 <- read.csv("data/2018_basil_deposition.csv", 
                    header = TRUE)

## QUICK VIZ ####
# Vetch 2018
boxplot((conspecific/(conspecific + other)) ~ treatment:site, 
        data = vetch18, 
        main = "Vetch 2018", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

# Vetch 2019
boxplot((conspecific/(conspecific + other)) ~ treatment:site, 
        data = vetch19, 
        main = "Vetch 2019", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

# Galium 2018
boxplot((conspecific/(conspecific + other)) ~ treatment:site, 
        data = galium18, 
        main = "Galium 2018", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

# Galium 2019
boxplot((conspecific/(conspecific + other)) ~ treatment:site, 
        data = galium19, 
        main = "Galium 2019", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

# Basil 2018
boxplot((conspecific/(conspecific + other)) ~ treatment, 
        data = basil18,
        main = "Basil 2018", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

## MODELS W/ BINOMIAL DISTRIBUTION ####
vet18_mod_bino <- glmmTMB(conspecific ~ treatment * site + (1|date),
                          data = vetch18,
                          family =  binomial(),
                          na.action = na.exclude, 
                          weights = (conspecific + other),
                          contrasts = list(treatment = "contr.sum",
                                           site = "contr.sum"))
summary(vet18_mod_bino)

## MODELS W/ BETA DISTRIBUTION ####
#### Vetch ####
# Keep only complete.cases
vetch18 <- vetch18[complete.cases(vetch18), ]

# Transform 1's in to .999
# for 2018
for(i in 1:nrow(vetch18)){
  if(vetch18$proportion_cons[i] == 1){
    vetch18$proportion_cons[i] <- 0.999
  }
}

# for 2019
for(i in 1:nrow(vetch19)){
  if(vetch19$proportion_cons[i] == 1){
    vetch19$proportion_cons[i] <- 0.999
  }
}

##### Models with Beta Distribution ####
###### 2018 Model Selection ####
vet18_mod <- glmmTMB(proportion_cons ~ treatment * site + (1|date),
                     data = vetch18,
                     family = beta_family(),
                     na.action = na.exclude, 
                     contrasts = list(treatment = "contr.sum",
                                      site = "contr.sum"))

# Diagnostics
summary(vet18_mod)
simulateResiduals(vet18_mod, # Good
                  plot = T)

vet18_mod2 <- glmmTMB(proportion_cons ~ treatment * site,
                      data = vetch18,
                      family = beta_family(),
                      na.action = na.exclude, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

# Compare
AIC(vet18_mod, 
    vet18_mod2) # without random effects better

# Diagnostics 
summary(vet18_mod2)
simulateResiduals(vet18_mod2, # Good
                  plot = T)

# LR Test
lrtest(vet18_mod2, 
       update(vet18_mod2, .~. -treatment:site)) # Good model

# Quick Viz
emmip(vet18_mod2, 
      treatment ~ site, 
      type = "response") # Effect driven by one plot difference


###### 2019 Model Selection ####
# Global Model
vet19_mod <- glmmTMB(proportion_cons ~ treatment * site + (1|date),
                     data = vetch19,
                     family = beta_family(),
                     na.action = na.exclude, 
                     contrasts = list(treatment = "contr.sum",
                                      site = "contr.sum"))

# Diagnostics
summary(vet19_mod)
simulateResiduals(vet19_mod, 
                  plot = T)

# Without random effect
vet19_mod2 <- glmmTMB(proportion_cons ~ treatment * site,
                      data = vetch19,
                      family = beta_family(),
                      na.action = na.exclude, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

# Compare
AIC(vet19_mod, 
    vet19_mod2) # simpler is better

# Diagnostics
summary(vet19_mod2)
simulateResiduals(vet19_mod2, # good
                  plot = T)

# Without interaction
vet19_mod3 <- glmmTMB(proportion_cons ~ treatment + site,
                      data = vetch19,
                      family = beta_family(),
                      na.action = na.exclude, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

# Compare
AIC(vet19_mod2, 
    vet19_mod3) # simpler is better

# Diagnostics
summary(vet19_mod3)
simulateResiduals(vet19_mod3, 
                  plot = T)

# LR Test
lrtest(vet19_mod3, 
       update(vet19_mod3, .~. -treatment))

# Quick Viz
emmip(vet19_mod3, 
      treatment ~ site, 
      type = "response")

### Galium ####
# keep only complete cases
galium18  <- galium18[complete.cases(galium18),]

# Transform 1's in to .999
# 2018
for(i in 1:nrow(galium18)){
  if(galium18$proportion_cons[i] == 1){
    galium18$proportion_cons[i] <- 0.999
  }
}

# 2019
for(i in 1:nrow(galium19)){
  if(galium19$proportion_cons[i] == 1){
    galium19$proportion_cons[i] <- 0.999
  }
}

##### Models with Beta Distribution ####
###### 2018 Model Selection #####
gal18_mod <- glmmTMB(proportion_cons ~ treatment * site + (1|date), 
                     data = galium18,
                     family = beta_family(),
                     na.action = na.exclude, 
                     contrasts = list(treatment = "contr.sum",
                                      site = "contr.sum"))

# Diagnostics
summary(gal18_mod)
simulateResiduals(gal18_mod, # Fails homogeneity test 
                  plot = T)
testZeroInflation(gal18_mod)

###### 2019 Model Selection #####
gal19_mod <- glmmTMB(proportion_cons ~ treatment * site + (1|date), 
                     data = galium19,
                     family = beta_family(),
                     na.action = na.exclude,
                     contrasts = list(treatment = "contr.sum",
                                      site = "contr.sum"))

# Diagnostics
summary(gal19_mod)
simulateResiduals(gal19_mod, 
                  plot = T)


#### Wild Basil ####
# NA's removed
basil18 <- basil18[complete.cases(basil18), ]

# Transform 1's to .999
for(i in 1:nrow(basil18)){
  if(basil18$proportion_cons[i] == 1){
    basil18$proportion_cons[i] <- 0.999
  }
}

##### 2018 Model Selection ####
bas18_mod <- glmmTMB(proportion_cons ~ treatment * subsample + (1|date), 
                     data = basil18,
                     family = beta_family(),
                     na.action = na.exclude, 
                     contrasts = list(treatment = "contr.sum"))

# Diagnostics
summary(bas18_mod)
simulateResiduals(bas18_mod, 
                  plot = T)