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