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


#### Vetch 2018 ####
# Quick Viz
boxplot((self/(self + other)) ~ treatment, 
        data = vetch18, 
        main = "Vetch", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

##### Models with beta distribution ####
# Keep only complete.cases
vetch18 <- vetch18[complete.cases(vetch18), ]

# Transform 1's in to .999
for(i in 1:nrow(vetch18)){
  if(vetch18$proportion_self[i] == 1){
    vetch18$proportion_self[i] <- 0.999
  }
}

# Global model vetch
vet_mod <- glmmTMB(proportion_self ~ treatment * site,
                   data = vetch18,
                   family = beta_family(),
                   na.action = na.exclude, 
                   contrasts = list(treatment = "contr.sum",
                                    site = "contr.sum"))

summary(vet_mod)
simulateResiduals(vet_mod, plot = T)
Anova(vet_mod)
lrtest(vet_mod, update(vet_mod, .~. -treatment:site))

emmip(vet_mod, treatment ~ site, type = "response")

#### Vetch 2019 ####
# Quick Viz
boxplot((self/(self + other)) ~ treatment, 
        data = vetch19, 
        main = "Vetch", 
        ylab = "Proportion of Conspecific Pollen", 
        ylim = range(0:1))

##### Models with beta distribution ####
# Keep only complete.cases
vetch18 <- vetch18[complete.cases(vetch18), ]

# Transform 1's in to .999
for(i in 1:nrow(vetch18)){
  if(vetch18$proportion_self[i] == 1){
    vetch18$proportion_self[i] <- 0.999
  }
}

# Global model vetch
vet_mod <- glmmTMB(proportion_self ~ treatment * site,
                   data = vetch18,
                   family = beta_family(),
                   na.action = na.exclude, 
                   contrasts = list(treatment = "contr.sum",
                                    site = "contr.sum"))

summary(vet_mod)
simulateResiduals(vet_mod, plot = T)
Anova(vet_mod)
lrtest(vet_mod, update(vet_mod, .~. -treatment:site))

emmip(vet_mod, treatment ~ site, type = "response")
