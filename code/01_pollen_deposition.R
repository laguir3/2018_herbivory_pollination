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
library(lme4)

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

vetch21 <- read.csv("data/2021_vetch_deposition.csv", 
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
ggplot(data = vetch18, 
       aes(y = proportion_cons,
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic() 
# NOTE: Deposition for the damage treatment, in the lower site, looks quite
# lower. Not sure the effect is there for the higher site. 

ggplot(data = vetch18, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(vetch18, 
                              subset = vetch18$treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = subset(vetch18, 
                               subset = vetch18$treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Overall, it looks like there are quite a few more stigmas with only 
# vetch pollen in the control plots. That's quite expected. 

# Vetch 2019
# boxplot
ggplot(data = vetch19, 
       aes(y = proportion_cons,
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic()
# NOTE: It looks like stigmas from damage plots received a lower proportion of 
# conspecific pollen in ALL SITES. Nice. 

# histogram
ggplot(data = vetch19, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(vetch19, 
                               subset = vetch19$treatment == "control"), 
                 fill = cb[4], 
                 alpha = 0.5) + 
  geom_histogram(data = subset(vetch19, 
                               subset = vetch19$treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) + 
  theme_classic()
# NOTE: This histogram also shows that stigmas in damage plots received lower 
# proportions of conspecific pollen. Much higher propability of 1's in the 
# control plots. 

# Vetch 2021
# Boxplot
ggplot(data = vetch21,
       aes(y = proportion_cons, 
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic()
# NOTE: This is unclear, there does seem to be a bit of a trend here as well,
# with lower conspecific deposition in the damage plots. However, I need to 
# add the rest of the data here. 

# Histogram
ggplot(data = vetch21, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(vetch21, 
                               subset = vetch21$treatment == "control"), 
                 fill = cb[4], 
                 alpha = 0.5) + 
  geom_histogram(data = subset(vetch21, 
                               subset = vetch21$treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) + 
  theme_classic()
# NOTE: This is much more unclear, possible effects a not clear at all.  

# Galium 2018
# Boxplot
ggplot(data = galium18,
       aes(y = proportion_cons, 
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic() + 
  ggtitle("2018 Galium Deposition")
# NOTE: There does see to be an effect at the lower site, with lower proportion 
# of conspecific pollen in teh damage plot. This pattern is not present in the 
# higher plot. 

# Histogram
ggplot(data = galium18, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(galium18, 
                               subset = galium18$treatment == "control"), 
                 fill = cb[4], 
                 alpha = 0.5) + 
  geom_histogram(data = subset(galium18, 
                               subset = galium18$treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Far more values close to one for the control stigmas. 

# Galium 2019
# Boxplot
ggplot(data = galium19,
       aes(y = proportion_cons, 
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic() + 
  ggtitle("2019 Galium Deposition")
# NOTE: This is a much more mixed batch of results. In gf1 and gf4, damage plots
# seem to have lower conspecific deposition. In gf2, control has lower 
# conspecific deposition, and in gf3 deposition seems to be about the same in 
# the control and damage plot. 

# Histogram 
ggplot(data = galium19, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(galium19, 
                               subset = galium19$treatment == "control"), 
                 fill = cb[4], 
                 alpha = 0.5) + 
  geom_histogram(data = subset(galium19, 
                               subset = galium19$treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Distributions look very, very similar. No clear effect. 

# Galium 2021
# NEED TO ENTER DATA

# Basil 2018
ggplot(data = basil18,
       aes(y = proportion_cons, 
           x = site, 
           color = treatment)) + 
  geom_boxplot() + 
  theme_classic() + 
  ggtitle("2018 Basil Deposition")
# NOTE: Definitely a difference, with a lower proportion of conspecific
# deposition in the damage plot. 

# Histogram
ggplot(data = basil18, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(basil18, 
                               subset = basil18$treatment == "control"), 
                 fill = cb[4], 
                 alpha = 0.5) + 
  geom_histogram(data = subset(basil18, 
                               subset = basil18$treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Clear difference, no 1's in the damage plots. Really great result. 


## MODELS W/ BETA DISTRIBUTION ####
#### Vetch ####
# Keep only complete.cases
vetch18 <- vetch18[complete.cases(vetch18), ]
vetch21 <- vetch21[complete.cases(vetch18), ]

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

for(i in 1:nrow(vetch21)){
  if(vetch21$proportion_cons[i] == 1){
    vetch21$proportion_cons[i] <- 0.999
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

## MODELS W/ ONE-INFLATED BETA DISTRIBUTION ####
## Vetch ####
vet18_oimodA <- glmmTMB(proportion_cons ~ treatment * site + (1|date),
                       data = subset(vetch18, 
                                     0 < proportion_cons & proportion_cons < 1),
                       family = beta_family(),
                       na.action = na.exclude, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))
vet18_oimodB <- glmmTMB((proportion_cons == 1) ~ treatment * site + (1|date),
                        data = vetch18,
                        family = binomial(),
                        na.action = na.exclude, 
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum"))

# Diagnostics
summary(vet18_oimodA)
summary(vet18_oimodB)
simulateResiduals(vet18_oimodA, # Good
                  plot = T)
simulateResiduals(vet18_oimodB, # Good
                  plot = T, as.factor = T)


vet18_oimodA <- glmmTMB(proportion_cons ~ treatment * site + (1|date),
                        data = subset(vetch18, 
                                      0 < proportion_cons & proportion_cons < 1),
                        family = binomial(),
                        na.action = na.exclude, 
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum"))



