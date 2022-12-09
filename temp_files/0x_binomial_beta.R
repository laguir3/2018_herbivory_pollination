#### DESCRIPTION ####
# This is the code for a question posted on StackOverflow regarding the
# interchangeability of an analysis with binomial and beta distributions for
# the pollen deposition. 
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



vetch18_full <- read.csv("data/2018_vetch_deposition.csv", 
                         header = TRUE)

vetch18_full <- vetch18_full[complete.cases(vetch18_full), ]

# vetch18 w/o 1 values
vetch18 <- vetch18_full %>% filter(proportion_cons < 1)

# vetch18 w/ 1 -> 0.999
vetch18_mod <- vetch18_full

for(i in 1:nrow(vetch18_mod)){
  if(vetch18_mod$proportion_cons[i] == 1){ # change 1's to 0.999
    vetch18_mod$proportion_cons[i] <- 0.999
  }
}


# histogram for full data
ggplot(data = vetch18_full, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(vetch18_full, 
                               subset = vetch18_full$treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = subset(vetch18_full, 
                               subset = vetch18_full$treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()


# histogram
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


# histogram w/ modified data
ggplot(data = vetch18_mod, 
       aes(x = proportion_cons)) +
  geom_histogram(data = subset(vetch18_mod, 
                               subset = vetch18_mod$treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = subset(vetch18_mod, 
                               subset = vetch18_mod$treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()

#### Models ####
# double column specification
bino_mod <- glmmTMB(cbind(conspecific, other) ~ treatment * site,
                    data = subset(vetch18_full),
                    family = binomial("logit"))
summary(bino_mod)
simulateResiduals(bino_mod, plot = T)


bino_mod2 <- glm(cbind(conspecific, other) ~ treatment * site,
                 data = subset(vetch18_full, 
                               proportion_cons < 1),
                 family = binomial("logit"))
summary(bino_mod2)
simulateResiduals(bino_mod2, plot = T)

# Beta model
beta_mod <- glmmTMB(proportion_cons ~ treatment * site,
                    data = subset(vetch18, 
                                  proportion_cons < 1),
                    family = beta_family())
summary(beta_mod)
simulateResiduals(beta_mod, plot = T)




