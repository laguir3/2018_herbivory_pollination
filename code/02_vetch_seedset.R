#### DESCRIPTION ####
# This is the analysis of the seedset data in Aguirre et al. Effects of
# herbivory on community pollination processes. 
# 
# Code by Luis A. Aguirre
# 
# #### LOAD ####
# Packages 
library(glmmTMB) 
library(tidyverse) # data wrangling
library(DHARMa) # model diagnostics
library(emmeans)
library(car) 
library(lmtest) # likelihood ratio test

# Set color-blind palette
cb <- c("#000000", # black
        "#E69F00", # orange
        "#56B4E9", # light blue
        "#009E73", # green
        "#F0E442", # yellow
        "#0072B2", # dark blue
        "#D55E00", # red
        "#CC79A7") # pink

# Load data
# 2018
seed18 <- read.csv("data/2018_vetch_seedset.csv",
                   header = T)

# 2021
seed21 <- read.csv("data/2021_vetch_seedset.csv", 
                   header = T)

#### QUICK VIZ ####
# 2018
ggplot(data = seed18, 
       aes(x = seeds)) + 
  geom_histogram(data = subset(seed18,
                               treatment == "control"),
                 fill = cb[4], 
                 alpha = 0.5) +
  geom_histogram(data = subset(seed18, 
                               treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) +
  theme_classic()
  
ggplot(data = seed21, 
       aes(x = seeds)) + 
  geom_histogram(data = subset(seed21,
                               treatment == "control"),
                 fill = cb[4], 
                 alpha = 0.5) +
  geom_histogram(data = subset(seed21, 
                               treatment == "damage"), 
                 fill = cb[7], 
                 alpha = 0.5) +
  theme_classic()

# TASKS
# change to 15 bins
# poisson models