#### DESCRIPTION ####
# Analysis for Aguirre et al. Effects of herbivory on community-wide 
# pollination processes.
# 
# Code by Luis A. Aguirre

### LOAD ####
# Packages 
library(glmmTMB) 
library(tidyverse)
library(DHARMa)
library(emmeans)
library(car)
library(lmtest)
library(lme4)
library(cowplot)
library(ggeffects)
library(patchwork)
library(flextable)
library(broom)
library(officer)

# Set color-blind palette
cb <- c("#000000", # black
        "#E69F00", # orange
        "#56B4E9", # light blue
        "#009E73", # green
        "#F0E442", # yellow
        "#0072B2", # dark blue
        "#D55E00", # red
        "#CC79A7") # pink

# Load common legend for all plots
common_legend <- readRDS("figures/common_legend.RDS")

# Function to change point shapes in plots created via ggeffects wrapper fun.
# Function to change point shapes in plots created via ggeffects wrapper fun.
manual_shape_change <- function(x, shape){
  # x = ggplot to modify
  temp <- ggplot_build(x)
  
  # v = vector of points to change
  v <- vector()
  if(nrow(temp$data[[1]]) == 4){
    v <- c(1,3)
  } else if(nrow(temp$data[[1]]) == 6) {
    v <- c(1,3,5)
  } else if(nrow(temp$data[[1]]) == 8){
    v <- c(1,3,5,7)
  }
  
  # set shape to triangle (17)
  temp$data[[1]]$shape[v] <- shape
  
  #
  return(as.ggplot(ggplot_gtable(temp)))
}

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

galium21 <- read.csv("data/2021_galium_deposition.csv", 
                     header = TRUE)

# Wild Basil
basil18 <- read.csv("data/2018_basil_deposition.csv", 
                    header = TRUE)

# Whorled Loosestrife
loose19 <- read.csv("data/2019_loosestrife_deposition.csv",
                    header = TRUE)

# Data format
# Treatment as factors
vetch18$treatment <- as.factor(vetch18$treatment)
vetch19$treatment <- as.factor(vetch19$treatment)
vetch21$treatment <- as.factor(vetch21$treatment)
galium18$treatment <- as.factor(galium18$treatment)
galium19$treatment <- as.factor(galium19$treatment)
galium21$treatment <- as.factor(galium21$treatment)
basil18$treatment <- as.factor(basil18$treatment)
loose19$treatment <- as.factor(loose19$treatment) 

# Site as factors
vetch18$plotpair <- as.factor(vetch18$site)
vetch19$plotpair <- as.factor(vetch19$site)
vetch21$plotpair <- as.factor(vetch21$site)
galium18$plotpair <- as.factor(galium18$site)
galium19$plotpair <- as.factor(galium19$site)
galium21$plotpair <- as.factor(galium21$site)
basil18$plotpair <- as.factor(basil18$site)
loose19$plotpair <- as.factor(loose19$site)

## QUICK VIZ ####
# Vetch 2018
ggplot(data = vetch18, 
       aes(y = proportion_het,
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2018 Vetch Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: Deposition for the damage treatment, in the lower site, looks quite
# lower. Not sure the effect is there for the higher site. 


ggplot(data = vetch18, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(vetch18, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(vetch18, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Overall, it looks like there are quite a few more stigmas with only 
# vetch pollen in the control plots. That's quite expected. 


# Vetch 2019
# boxplot
ggplot(data = vetch19, 
       aes(y = proportion_het,
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2019 Vetch Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: It looks like stigmas from damage plots received a lower proportion of 
# conspecific pollen in ALL SITES. Nice. 

# histogram
ggplot(data = vetch19, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(vetch19, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(vetch19, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: This histogram also shows that stigmas in damage plots received lower 
# proportions of conspecific pollen. Much higher propability of 1's in the 
# control plots. 

# Vetch 2021
# Boxplot
ggplot(data = vetch21,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() +
  ggtitle("2021 Vetch Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: This is unclear, there does seem to be a bit of a trend here as well,
# with lower conspecific deposition in the damage plots. However, I need to 
# add the rest of the data here. 

# Histogram
ggplot(data = vetch21, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(vetch21, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(vetch21, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: This is much more unclear, possible effects a not clear at all.  

# Galium 2018
# Boxplot
ggplot(data = galium18,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2018 Galium Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: There does see to be an effect at the lower site, with lower proportion 
# of conspecific pollen in teh damage plot. This pattern is not present in the 
# higher plot. 

# Histogram
ggplot(data = galium18, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(galium18, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(galium18, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Far more values close to one for the control stigmas. 

# Galium 2019
# Boxplot
ggplot(data = galium19,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2019 Galium Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: This is a much more mixed batch of results. In gf1 and gf4, damage plots
# seem to have lower conspecific deposition. In gf2, control has lower 
# conspecific deposition, and in gf3 deposition seems to be about the same in 
# the control and damage plot. 

# Histogram 
ggplot(data = galium19, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(galium19, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(galium19, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Distributions look very, very similar. No clear effect. 

# Galium 2021
# Boxplot
ggplot(data = galium21,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2021 Galium Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")

# Histogram 
ggplot(data = galium21, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(galium21, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(galium21, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Distributions look very, very similar. No clear effect. 


# Basil 2018
ggplot(data = basil18,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2018 Basil Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: Definitely a difference, with a lower proportion of conspecific
# deposition in the damage plot. 

# Histogram
ggplot(data = basil18, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(basil18, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(basil18, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Clear difference, no 1's in the damage plots. Really great result. 

# Whorled Loosestrife 2018
ggplot(data = loose19,
       aes(y = proportion_het, 
           x = plotpair, 
           fill = treatment)) + 
  geom_violin() + 
  theme_classic() + 
  ggtitle("2019 Loosestrife Pollen Deposition") +
  ylab("Proportion of Conspecific Pollen")
# NOTE: Definitely a difference, with a lower proportion of conspecific
# deposition in the damage plot. 

# Histogram
ggplot(data = loose19, 
       aes(x = proportion_het)) +
  geom_histogram(data = filter(loose19, 
                               treatment == "control"), 
                 fill = cb[4], # green
                 alpha = 0.5) + 
  geom_histogram(data = filter(loose19, 
                               treatment == "damage"), 
                 fill = cb[7], # red
                 alpha = 0.5) + 
  theme_classic()
# NOTE: Clear difference, no 1's in the damage plots. Really great result.

## MODELS W/ BETA DISTRIBUTION ####
## Create list for table summaries
app_summs <- vector(mode = "list", 
                    length = 8) # number of pollen dep. models

## Create list for Likelihood Ratio Test (LRT's)
depo_lrts <- vector(mode = "list", 
               length = 8) # number of pollen dep. models

# name list objects
names(app_summs) <- c("2018 Vetch Deposition", 
                      "2019 Vetch Deposition", 
                      "2021 Vetch Deposition", 
                      "2018 Galium Deposition", 
                      "2019 Galium Deposition", 
                      "2021 Galium Deposition", 
                      "2018 Wild Basil Deposition", 
                      "2019 Whorled Loosestrife Deposition")

names(depo_lrts) <- names(app_summs)


### Vetch ####
# Keep only complete.cases
vetch18 <- vetch18[complete.cases(vetch18), ]
vetch21 <- vetch21[complete.cases(vetch18), ]

##### 2018 Model Selection ####
vet18_mod <- glmmTMB(proportion_het ~ treatment * plotpair +
                       (1|date),
                     data = vetch18,
                     family = beta_family(),
                     na.action = na.exclude, 
                     ziformula = ~ treatment * plotpair,
                     contrasts = list(treatment = "contr.sum",
                                      plotpair = "contr.sum"))

# Model check
summary(vet18_mod)
simulateResiduals(vet18_mod, # Good
                  plot = T)

# Model w/o random effects
vet18_mod2 <- glmmTMB(proportion_het ~ treatment * plotpair,
                      data = vetch18,
                      family = beta_family(),
                      na.action = na.exclude, 
                      ziformula = ~ treatment * plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check 
summary(vet18_mod2)
simulateResiduals(vet18_mod2, # Good
                  plot = T)

# Compare
AIC(vet18_mod, 
    vet18_mod2) # Simpler is better

# Model w/o interaction
vet18_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair,
                      data = vetch18,
                      family = beta_family(),
                      na.action = na.exclude, 
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Compare
AIC(vet18_mod, 
    vet18_mod2, # Model 2 is best
    vet18_mod3)

# Model check
summary(vet18_mod3)
simulateResiduals(vet18_mod3, # Good
                  plot = T)


###### Best Fit ####
app_summs[[1]] <- vet18_mod2

# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment * plotpair
# Zero inflation:                  ~treatment * plotpair
# Data: vetch18
# 
# AIC      BIC   logLik deviance df.resid 
# -74.5    -47.0     46.3    -92.5      148 
# 
# 
# Dispersion parameter for beta family (): 8.53 
# 
# Conditional model:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -1.732445   0.082940 -20.888  < 2e-16 ***
# treatment1           -0.006775   0.074029  -0.092 0.927083    
# plotpair1            -0.020404   0.074032  -0.276 0.782850    
# treatment1:plotpair1 -0.246177   0.074447  -3.307 0.000944 ***
# 
# Zero-inflation model:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -1.6868     0.3011  -5.602 2.11e-08 ***
# treatment1         0.9203     0.3011   3.057  0.00224 ** 
# plotpair1         -0.3264     0.3011  -1.084  0.27825    
# treatment1:plotpair1   0.7300     0.3011   2.425  0.01532 *  

###### LRT ####
depo_lrts[[1]] <- lrtest(vet18_mod2, 
                    update(vet18_mod2, 
                           .~. -treatment:plotpair)) # Good model
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment * plotpair
# Model 2: proportion_het ~ treatment + plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)   
# 1   9 46.273                        
# 2   8 40.880 -1 10.786   0.001023 **



# Dertermine difference between treatments
pairs(emmeans(vet18_mod2, 
              ~ treatment * plotpair, 
              type = "response"))
# contrast                  odds.ratio    SE  df null t.ratio p.value
# control gf1 / damage gf1       0.603 0.130 148    1  -2.346  0.0925
#  ***
# control gf3 / damage gf3       1.614 0.329 148    1   2.346  0.0925
# NOTE: Inflated because these are tukey tests

print(ggemmeans(vet18_mod, 
                ~ treatment * plotpair, 
                type = "fe.zi"), 
      digits = 5)
# treatment plotpair response     SE  df lower.CL upper.CL
# control   gf1     0.119 0.0191 147   0.0858    0.162
# damage    gf1     0.182 0.0189 147   0.1480    0.223
# control   gf3     0.187 0.0220 147   0.1469    0.234
# damage    gf3     0.124 0.0166 147   0.0951    0.161
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
vet18_plot <- plot(ggpredict(vet18_mod2, 
                             ~ plotpair + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
vet18_plot <- vet18_plot + 
  ggtitle(expression(paste("A) 2018 - ", 
                           italic("V. cracca")))) + 
  labs(# Turn on/off
    # x = "Plot Pair",  
    x = "",
    y = "Proportion of \n Heterospecific Pollen"
  ) +
  ylim(range(0,.30)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) +
  scale_x_discrete(limits=c("A", "C")) +
  theme(legend.position = "none") +
  annotate(geom = "text",
           x = 1.6,
           y = .27,
           label = paste("Treatment: P = 0.93\n",
                         "Treatment x Plot Pair: \n P = 0.01"), 
           size = 8) 

# change shape for treatment
vet18_plot <- manual_shape_change(vet18_plot, 17)


##### 2019 Model Selection ####

# Global Model
vet19_mod <- glmmTMB(proportion_het ~ treatment + plotpair + vetch +
                       (1|date),
                     data = vetch19[complete.cases(vetch19),], # ignore na's
                     family = beta_family(),
                     ziformula = ~ treatment + plotpair,
                     contrasts = list(treatment = "contr.sum",
                                      plotpair = "contr.sum"))

# Model check
summary(vet19_mod)
simulateResiduals(vet19_mod, 
                  plot = T) # A lot of dispersion and fails KS test
testDispersion(vet19_mod)

# Model w/ random effects
vet19_mod2 <- glmmTMB(proportion_het ~ treatment + plotpair + vetch,
                      data = vetch19[complete.cases(vetch19),], # ignore na's
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Compare
AIC(vet19_mod, 
    vet19_mod2) # simpler is better

# Model check
summary(vet19_mod2)
simulateResiduals(vet19_mod2,
                  plot = T) # Still lots of dispersion
testDispersion(vet19_mod2)

# Without resources
vet19_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair,
                      data = vetch19[complete.cases(vetch19),], #ignore na's
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Compare
AIC(vet19_mod, 
    vet19_mod2, 
    vet19_mod3) # better, but overdispersion possible

# Model check
summary(vet19_mod3)
simulateResiduals(vet19_mod3, 
                  plot = T) # Seems like this is overdispersed


# Use complete dataset
refit_vet19_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair,
                            data = vetch19, # all observations
                            family = beta_family(),
                            ziformula = ~ treatment + plotpair,
                            contrasts = list(treatment = "contr.sum",
                                             plotpair = "contr.sum"))

# Model check
summary(refit_vet19_mod3)
simulateResiduals(refit_vet19_mod3, 
                  plot = T) # Seems like this is overdispersed

###### Best Fit ####
app_summs[[2]] <- refit_vet19_mod3
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment + plotpair
# Zero inflation:                  ~treatment + plotpair
# Data: vetch19
# 
# AIC      BIC   logLik deviance df.resid 
# -288.5   -248.9    155.3   -310.5      261 
# 
# 
# Dispersion parameter for beta family (): 14.4 
# 
# Conditional model:
#            Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.22862    0.06267  -35.56  < 2e-16 ***
# treatment1  -0.16849    0.05227   -3.22  0.00127 ** 
# plotpair1        0.23941    0.08343    2.87  0.00411 ** 
# plotpair2       -0.01324    0.10564   -0.13  0.90024    
# plotpair3       -0.16057    0.08493   -1.89  0.05866 .  
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5957     0.1985  -8.039 9.05e-16 ***
# treatment1    0.5078     0.1625   3.124 0.001781 ** 
# plotpair1        -0.1581     0.2839  -0.557 0.577610    
# plotpair2        -0.4455     0.3908  -1.140 0.254334    
# plotpair3         0.8318     0.2383   3.491 0.000481 ***


###### LRT ####
depo_lrts[[2]] <- lrtest(vet19_mod3, 
                    update(vet19_mod3, 
                           .~. -treatment))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + plotpair
# Model 2: proportion_het ~ plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)   
# 1  11 155.26                        
# 2  10 150.17 -1 10.186   0.001415 **

# Determine pairwise differences
pairs(emmeans(vet19_mod3, 
              "treatment", 
              type = "response"))
# contrast         odds.ratio     SE  df null t.ratio p.value
# control / damage        0.7 0.0781 241    1  -3.193  0.0016

emmeans(vet19_mod3, 
        ~ treatment + plotpair, 
        type = "response")
# treatment plotpair response      SE  df lower.CL upper.CL
# control   gf1    0.1036 0.01132 241   0.0834   0.1281
# damage    gf1    0.1417 0.01481 241   0.1150   0.1735
# control   gf2    0.0823 0.01135 241   0.0625   0.1076
# damage    gf2    0.1135 0.01382 241   0.0890   0.1437
# control   gf3    0.0721 0.00770 241   0.0583   0.0888
# damage    gf3    0.0999 0.00994 241   0.0819   0.1212
# control   gf4    0.0786 0.01087 241   0.0596   0.1028
# damage    gf4    0.1085 0.01389 241   0.0841   0.1391
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale

###### Plot ####
vet19_plot <- plot(ggpredict(vet19_mod3, 
                             ~ plotpair + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
vet19_plot <- vet19_plot + 
  ggtitle(expression(paste("B) 2019 - ", 
                           italic("V. cracca")))) +
  labs(x = "Plot Pair",
       # Turn on/off
       # y = "Proportion of \n Heterospecific Pollen", 
       y = ""
  ) +
  ylim(range(0,.30)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("A", 
                              "B",
                              "C", 
                              "D")) +
  theme(legend.position = "none") +
  annotate(geom = "text", 
           y = .3, 
           x = 3, 
           label = "Treatment: P = 0.01", 
           size = 8)

# change shape for treatment
vet19_plot <- manual_shape_change(vet19_plot, 17)


##### 2021 Model Selection ####
# Global Model with vetch as resources
vet21_modv <- glmmTMB(proportion_het ~ treatment * plotpair + vetch +  
                        (1|date),
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(vet21_modv)
simulateResiduals(vet21_modv,
                  plot = T)

# Global Model with milkweed as resources
vet21_modm <- glmmTMB(proportion_het ~ treatment * plotpair + milk_in +  
                        (1|date),
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(vet21_modm)
simulateResiduals(vet21_modm,
                  plot = T)

# Compare
AIC(vet21_modm, # Essentially the same
    vet21_modv)

# Without resources
vet21_mod2 <- glmmTMB(proportion_het ~ treatment * plotpair +  
                        (1|date),
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(vet21_mod2)
simulateResiduals(vet21_mod2,
                  plot = T)

# Compare
AIC(vet21_modv, 
    vet21_mod2) # Better

# Model w/o random effects
vet21_mod3 <- glmmTMB(proportion_het ~ treatment * plotpair ,
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))



# Model check
summary(vet21_mod3)
simulateResiduals(vet21_mod3,
                  plot = T)

# Compare
AIC(vet21_mod2, 
    vet21_mod3) # simples is better

# Model w/o interaction
vet21_mod4 <- glmmTMB(proportion_het ~ treatment + plotpair,
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(vet21_mod4)
simulateResiduals(vet21_mod4,
                  plot = T)

#Compare
AIC(vet21_mod3,
    vet21_mod4) # Best


# Model w/o plotpair
vet21_mod5 <- glmmTMB(proportion_het ~ treatment,
                      data = vetch21,
                      family = beta_family(),
                      ziformula = ~ treatment,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(vet21_mod5)
simulateResiduals(vet21_mod5,
                  plot = T)

# Model modifying to keep only sig terms in the fe and zi components
vet21_mod5b <- glmmTMB(proportion_het ~ treatment + plotpair,
                       data = vetch21,
                       family = beta_family(),
                       ziformula = ~ plotpair,
                       contrasts = list(treatment = "contr.sum",
                                        plotpair = "contr.sum"))

# Model check
summary(vet21_mod5b)
simulateResiduals(vet21_mod5b,
                  plot = T)

# Compare
AIC(vet21_mod4, 
    vet21_mod5, 
    vet21_mod5b) # Best

###### Best Fit ####
app_summs[[3]] <- vet21_mod5b
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment + plotpair
# Zero inflation:                  ~plotpair
# Data: vetch21
# 
# AIC      BIC   logLik deviance df.resid 
# -153.4   -127.5     86.7   -173.4       89 
# 
# 
# Dispersion parameter for beta family ():   59 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -3.63996    0.11885 -30.626   <2e-16 ***
#   treatment1  -0.19663    0.09605  -2.047   0.0406 *  
#   plotpair1       -0.23039    0.20295  -1.135   0.2563    
# plotpair2        0.29348    0.19422   1.511   0.1308    
# plotpair3        0.03615    0.14416   0.251   0.8020    
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.2324     0.2292  -1.014   0.3106  
# plotpair1         0.4330     0.3918   1.105   0.2691  
# plotpair2         0.8514     0.4030   2.113   0.0346 *
#   plotpair3        -0.1305     0.3248  -0.402   0.6878  

###### LRT #####
depo_lrts[[3]] <- lrtest(vet21_mod5b, 
                    update(vet21_mod5b,
                           .~. -treatment))

# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment
# Model 2: proportion_het ~ 1
# #Df LogLik Df  Chisq Pr(>Chisq)  
# 1   5 80.643                       
# 2   4 78.961 -1 3.3652    0.06659 . # Fail to reject null


# Dertermine difference between treatments
pairs(emmeans(vet21_mod5b, 
              "treatment", 
              type = "response"))
# contrast         odds.ratio    SE df null t.ratio p.value
# control / damage      0.699 0.136 92    1  -1.844  0.0684


print(ggpredict(vet21_mod5b, 
                terms = c("treatment", "plotpair"), 
                type = "fe"), 
      digits = 5)
# treatment plotpair response      SE df lower.CL upper.CL
# control   gf1    0.0168 0.00476 88  0.00958   0.0295
# damage    gf1    0.0248 0.00634 88  0.01484   0.0410
# control   gf2    0.0281 0.00698 88  0.01711   0.0458
# damage    gf2    0.0411 0.00991 88  0.02534   0.0660
# control   gf3    0.0219 0.00418 88  0.01493   0.0319
# damage    gf3    0.0321 0.00527 88  0.02310   0.0444
# control   gf4    0.0192 0.00427 88  0.01228   0.0298
# damage    gf4    0.0281 0.00549 88  0.01904   0.0414
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
vet21_plot <- plot(ggpredict(vet21_mod5b, 
                             ~ plotpair + treatment, 
                             type = "fe"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
vet21_plot <- vet21_plot + 
  ggtitle(expression(paste("C) 2021 - ", 
                           italic("V. cracca")))) +
  labs(# Turn on/off
    # x = "Plot Pair",  
    x = "",
    # y = "Proportion of \n Heterospecific Pollen", 
    y = ""
  ) +
  ylim(range(0,.10)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("A", 
                              "B",
                              "C", 
                              "D")) +
  theme(legend.position = "none", 
        plot.title.position = ) +
  annotate(geom = "text", 
           y = .1, 
           x = 3,
           label = "Treatment: P = 0.04", 
           size = 8) 

# change shape for treatment
vet21_plot <- manual_shape_change(vet21_plot, 17)

##### All Vetch Plots ####
vetch_plots <- (vet18_plot + vet19_plot + vet21_plot) / 
  common_legend + 
  plot_annotation(expression(paste(italic("Vicia cracca "), 
                                   "Pollen Deposition"))) + 
  plot_layout(guides = "collect", 
              height = c(8,1)) &
  theme(plot.title = element_text(size = 30)) 

# Save
ggsave("figures/vetch_deposition.png", 
       last_plot(), 
       device = "png",
       width = 18, 
       height = 7, 
       units = "in",
       dpi = 300)

### Galium ####
# keep only complete cases
galium18  <- galium18[complete.cases(galium18),]


##### 2018 Model Selection #####
# Global Model
gal18_mod <- glmmTMB(proportion_het ~ treatment * plotpair +
                       (1|date), 
                     data = galium18,
                     family = beta_family(),
                     ziformula = ~ treatment * plotpair,
                     contrasts = list(treatment = "contr.sum",
                                      plotpair = "contr.sum"))

# Model check
summary(gal18_mod)
simulateResiduals(gal18_mod, 
                  plot = T)

# w/o random effects
gal18_mod2 <- glmmTMB(proportion_het ~ treatment * plotpair, 
                      data = galium18,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal18_mod2)
simulateResiduals(gal18_mod2, 
                  plot = T)

# Compare 
AIC(gal18_mod, 
    gal18_mod2) # Better w/o random effects

# w/o interactions
gal18_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair, 
                      data = galium18,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Compare
AIC(gal18_mod, 
    gal18_mod2, 
    gal18_mod3)# Best

# Model check
summary(gal18_mod3)
simulateResiduals(gal18_mod3, 
                  plot = T)

# Without Treatment
gal18_mod4 <- glmmTMB(proportion_het ~ plotpair, 
                      data = galium18,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal18_mod4)
simulateResiduals(gal18_mod4, 
                  plot = T)

AIC(gal18_mod, 
    gal18_mod2, 
    gal18_mod3, 
    gal18_mod4) # Best


###### Best Fit ####
app_summs[[4]] <- gal18_mod3
# Family: beta  ( logit )
# Formula:          proportion_het ~ plotpair
# Zero inflation:                  ~treatment * plotpair
# Data: galium18
# 
# AIC      BIC   logLik deviance df.resid 
# -11.8      9.6     12.9    -25.8      152 
# 
# 
# Dispersion parameter for beta family (): 13.2 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.0831     0.0917  -22.72   <2e-16 ***
#   plotpair1        -0.0361     0.0814   -0.44     0.66    
# 
# Zero-inflation model:
#                 Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -0.134      0.180   -0.75  0.45531    
# treatment1          0.762      0.180    4.23  2.3e-05 ***
# plotpair1               0.118      0.180    0.66  0.51031    
# treatment1:plotpair1    0.609      0.180    3.38  0.00071 ***

###### LRT ####
depo_lrts[[4]] <- lrtest(gal18_mod3, 
                    update(gal18_mod2, 
                           .~. -treatment)) # without interaction

# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment * plotpair
# Model 2: proportion_het ~ treatment + plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)
# 1   9 13.202                     
# 2   8 12.920 -1 0.5644     0.4525 # No main effect of interaction, but it 
#                                   # matters for zero inflation

# Dertermine difference between treatments
pairs(emmeans(gal18_mod3, 
              "treatment", 
              type = "response"))
# contrast         odds.ratio    SE  df null t.ratio p.value
# control / damage      0.974 0.183 150    1  -0.140  0.8885
# 
# Results are averaged over the levels of: plotpair 
# Tests are performed on the log odds ratio scale 

print(ggemmeans(gal18_mod3, 
                ~ treatment + plotpair, 
                type = "re.zi"), 
      digits = 3)
# treatment plotpair response     SE  df lower.CL upper.CL
# control   gf1    0.0951 0.0240 150   0.0571    0.154
# damage    gf1    0.1103 0.0136 150   0.0862    0.140
# control   gf3    0.1205 0.0174 150   0.0901    0.159
# damage    gf3    0.1090 0.0154 150   0.0820    0.143
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
gal18_plot <- plot(ggpredict(gal18_mod3, 
                             ~ plotpair + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
gal18_plot <- gal18_plot + 
  ggtitle(expression(paste("D) 2018 - ", 
                           italic("G. palustre")))) +
  labs(# Turn on/off
    # x = "Plot Pair",  
    x = "",
    y = "Proportion of \n Heterospecific Pollen"
  ) +
  ylim(range(0,.15)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits=c("A", "C")) +
  theme(legend.position = "none") +
  annotate(geom = "text",
           y = .15, 
           x = 1.75, 
           label = "Treatment: P = 0.94", 
           size = 8)

# change shape for treatment
gal18_plot <- manual_shape_change(gal18_plot, 17)


##### 2019 Model Selection #####
# Global model with vetch as resource
gal19_modv <- glmmTMB(proportion_het ~ treatment * plotpair + vetch +
                        (1|date), 
                      data = galium19[complete.cases(galium19),], # ignore na's,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_modv)
simulateResiduals(gal19_modv, 
                  plot = T)

# Global model with milkweed as resource
gal19_modm <- glmmTMB(proportion_het ~ treatment * plotpair + milk_in +
                        (1|date), 
                      data = galium19[complete.cases(galium19),], # ignore na's,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_modm)
simulateResiduals(gal19_modm, 
                  plot = T)

# Global model with vetch as resource
gal19_modg <- glmmTMB(proportion_het ~ treatment * plotpair + galium +
                        (1|date), 
                      data = galium19[complete.cases(galium19),], # ignore na's,
                      family = beta_family(),
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_modg)
simulateResiduals(gal19_modg, 
                  plot = T)

# Compare
AIC(gal19_modv,
    gal19_modm, 
    gal19_modg) # Ther are essentially all the same

# Without date
gal19_mod2 <- glmmTMB(proportion_het ~ treatment * plotpair + galium, 
                      data = galium19[complete.cases(galium19),], # ignore na's
                      family = beta_family(), 
                      ziformula = ~ treatment * plotpair,
                      contrasts = list(treatment = "contr.sum", 
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_mod2)
simulateResiduals(gal19_mod2, 
                  plot = T)

# Compare
AIC(gal19_modm,
    gal19_mod2) # better

# Without interactions
gal19_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair + galium, 
                      data = galium19[complete.cases(galium19),], # ignore na's
                      family = beta_family(), 
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum", 
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_mod3)
simulateResiduals(gal19_mod3, 
                  plot = T)

# Compare
AIC(gal19_mod2, 
    gal19_mod3) # simpler is better

# Without resources
gal19_mod4 <- glmmTMB(proportion_het ~ treatment + plotpair, 
                      data = galium19[complete.cases(galium19),], # ignore na's
                      family = beta_family(), 
                      ziformula = ~ treatment + plotpair,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_mod4)
simulateResiduals(gal19_mod4,
                  plot = T)

# Compare
AIC(gal19_mod3, 
    gal19_mod4) # simpler is better

# refit with full data set
refit_gal19_mod4 <- glmmTMB(proportion_het ~ treatment + plotpair, 
                            data = galium19,
                            family = beta_family(), 
                            ziformula = ~ treatment + plotpair,
                            contrasts = list(treatment = "contr.sum",
                                             plotpair = "contr.sum"))

# Model check
summary(refit_gal19_mod4)
simulateResiduals(refit_gal19_mod4,
                  plot = T) # some outlier, but summaries are close enough

# Without zero inflation
gal19_mod5 <- glmmTMB(proportion_het ~ treatment + plotpair, 
                      data = galium19,
                      family = beta_family(), 
                      ziformula = ~ 1,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_mod5)
simulateResiduals(gal19_mod5,
                  plot = T)

# Compare 
AIC(refit_gal19_mod4, 
    gal19_mod5)

# Without Treatment
gal19_mod6 <- glmmTMB(proportion_het ~  plotpair, 
                      data = galium19,
                      family = beta_family(), 
                      ziformula = ~ 1,
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal19_mod6)
simulateResiduals(gal19_mod6,
                  plot = T)

# Compare 
AIC(gal19_mod5,
    gal19_mod6) # Simpler is better

###### Best Fit #####
app_summs[[5]] <- refit_gal19_mod4
# Family: beta  ( logit )
# Formula:          proportion_het ~ plotpair
# Zero inflation:                  ~1
# Data: galium19
# 
# AIC      BIC   logLik deviance df.resid 
# -336.0   -315.3    174.0   -348.0      226 
# 
# 
# Dispersion parameter for beta family (): 18.7 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.59817    0.06853  -37.91  < 2e-16 ***
#   plotpair1       -0.25942    0.10086   -2.57  0.01011 *  
#   plotpair2        0.26467    0.09229    2.87  0.00413 ** 
#   plotpair3        0.07037    0.09496    0.74  0.45869    
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1688     0.1544  -7.571  3.7e-14 ***

###### LRT ####
depo_lrts[[5]] <- lrtest(gal19_mod5, 
                    update(gal19_mod5, 
                           .~. -treatment))

# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + plotpair
# Model 2: proportion_het ~ plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)
# 1   7  174.2                     
# 2   6  174.0 -1 0.3945       0.53

# Determine effect of damage
pairs(emmeans(gal19_mod5, 
              "treatment", 
              type = "response"))
# contrast         odds.ratio    SE  df null t.ratio p.value
# control / damage       1.07 0.122 225    1   0.628  0.5305
# 
# Results are averaged over the levels of: plotpair 
# Tests are performed on the log odds ratio scale 

emmeans(gal19_mod5, 
        ~ treatment + plotpair, 
        type = "response")
# treatment plotpair response      SE  df lower.CL upper.CL
# control   gf1    0.0559 0.00698 225   0.0436   0.0714
# damage    gf1    0.0523 0.00688 225   0.0403   0.0676
# control   gf2    0.0913 0.00976 225   0.0738   0.1124
# damage    gf2    0.0856 0.00925 225   0.0690   0.1056
# control   gf3    0.0764 0.00870 225   0.0609   0.0954
# damage    gf3    0.0715 0.00824 225   0.0569   0.0896
# control   gf4    0.0666 0.00992 225   0.0495   0.0890
# damage    gf4    0.0623 0.00948 225   0.0461   0.0838
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
gal19_plot <- plot(ggpredict(gal19_mod5, 
                             ~ plotpair + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
gal19_plot <- gal19_plot + 
  ggtitle(expression(paste("E) 2019 - ", 
                           italic("G. palustre")))) +
  labs(# Turn on/off
    x = "Plot Pair",
    # y = "Proportion of \n Heterospecific Pollen", 
    y = ""
  ) +
  ylim(range(0,.15)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("A", 
                              "B",
                              "C", 
                              "D")) +
  theme(legend.position = "none") +
  annotate(geom = "text",
           y = .15, 
           x = 3, 
           label = "Treatment: P = 0.53", 
           size = 8)

# change shape for treatment
gal19_plot <- manual_shape_change(gal19_plot, 17)


##### 2021 Model Selection ####
# Global Model with vetch as resource
gal21_modv <- glmmTMB(proportion_het ~ treatment + plotpair + vetch + 
                        (1|date), 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))
# NOTE: won't converge with interactions.

# Model check
summary(gal21_modv)
simulateResiduals(gal21_modv, 
                  plot = T)

# Global Model with milkweed as resource
gal21_modm <- glmmTMB(proportion_het ~ treatment + plotpair + milk_in +
                        (1|date), 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))
# NOTE: won't converge with interactions.

# Model check
summary(gal21_modm)
simulateResiduals(gal21_modm, 
                  plot = T)

# Global Model with galium as resource
gal21_modg <- glmmTMB(proportion_het ~ treatment + plotpair + galium + 
                        (1|date), 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))
# NOTE: won't converge with interactions. Zero inflation part looks funky

# Model check
summary(gal21_modg)
simulateResiduals(gal21_modg, 
                  plot = T)

# Compare
AIC(gal21_modv, 
    gal21_modm, 
    gal21_modg) # all the same

# Without random effects
gal21_mod2 <- glmmTMB(proportion_het ~ treatment + plotpair + galium, 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal21_mod2)
simulateResiduals(gal21_mod2, 
                  plot = T)

# Compare
AIC(gal21_modg,
    gal21_mod2) # Best

# Without Resource
gal21_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair, 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ treatment + plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal21_mod3)
simulateResiduals(gal21_mod3, 
                  plot = T)

# Compare
AIC(gal21_mod2, 
    gal21_mod3) # simpler

# Without Treatment
gal21_mod4 <- glmmTMB(proportion_het ~ plotpair, 
                      data = galium21,
                      family = beta_family(),
                      ziformula = ~ plotpair, 
                      contrasts = list(treatment = "contr.sum",
                                       plotpair = "contr.sum"))

# Model check
summary(gal21_mod4)
simulateResiduals(gal21_mod4, 
                  plot = T)

# Compare
AIC(gal21_mod3, 
    gal21_mod4) # better

###### Best Fit ####
app_summs[[6]] <- gal21_mod3
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment + plotpair
# Zero inflation:                  ~ treatment + plotpair
# Data: galium21
# 
# AIC      BIC   logLik deviance df.resid 
# 40.3     68.6     -9.1     18.3       86 
# 
# 
# Dispersion parameter for beta family (): 16.7 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.6115     0.2726   -9.58   <2e-16 ***
# treatment1   -0.0868     0.1514   -0.57     0.57    
# plotpair1        -0.6892     0.6860   -1.00     0.32    
# plotpair2         0.2092     0.3085    0.68     0.50    
# plotpair3         0.1002     0.3024    0.33     0.74    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.1907     0.3232    3.68  0.00023 ***
# treatment1    0.1674     0.2455    0.68  0.49524    
# plotpair1         1.8255     0.7944    2.30  0.02157 *  
# plotpair2        -1.3192     0.4734   -2.79  0.00532 ** 
# plotpair3        -0.0851     0.4139   -0.21  0.83711    

###### LRT ####
depo_lrts[[6]] <- lrtest(gal21_mod3, 
                    update(gal21_mod3, 
                           .~. -treatment))

# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + plotpair
# Model 2: proportion_het ~ plotpair
#   #Df  LogLik Df  Chisq Pr(>Chisq)
# 1  11 -9.1467                     
# 2  10 -9.3108 -1 0.3282     0.5667 # fail to reject null

# Determine effect of damage
pairs(emmeans(gal21_mod3, 
              "treatment", 
              type = "response"))
# contrast         odds.ratio    SE  df null t.ratio p.value
# control / damage       1.07 0.122 225    1   0.628  0.5305
# 
# Results are averaged over the levels of: plotpair 
# Tests are performed on the log odds ratio scale 

emmeans(gal21_mod3, 
        ~ treatment + plotpair, 
        type = "response")
# treatment plotpair response      SE  df lower.CL upper.CL
# control   gf1    0.0559 0.00698 225   0.0436   0.0714
# damage    gf1    0.0523 0.00688 225   0.0403   0.0676
# control   gf2    0.0913 0.00976 225   0.0738   0.1124
# damage    gf2    0.0856 0.00925 225   0.0690   0.1056
# control   gf3    0.0764 0.00870 225   0.0609   0.0954
# damage    gf3    0.0715 0.00824 225   0.0569   0.0896
# control   gf4    0.0666 0.00992 225   0.0495   0.0890
# damage    gf4    0.0623 0.00948 225   0.0461   0.0838
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
gal21_plot <- plot(ggpredict(gal21_mod3, 
                             ~ plotpair + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
gal21_plot <- gal21_plot + 
  ggtitle(expression(paste("F) 2021 - ", 
                           italic("G. palustre")))) +
  labs(# Turn on/off
    # x = "Plot Pair",  
    x = "",
    # y = "Proportion of \n Heterospecific Pollen", 
    y = ""
  ) +
  ylim(range(0,.15)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("A", 
                              "B",
                              "C", 
                              "D")) +
  theme(legend.position = "none") +
  annotate(geom = "text", 
           y = .15, 
           x = 3, 
           label = "Treatment: P = 0.57", 
           size = 8)

# change shape for treatment
gal21_plot <- manual_shape_change(gal21_plot, 17)


##### All Galium Plots ####
galium_plots <- (gal18_plot + gal19_plot + gal21_plot) /
  common_legend +
  plot_annotation(expression(paste(italic("Galium palustre "),
                                   "Pollen Deposition"))) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 30))


# Save
ggsave("figures/galium_deposition.png", 
       last_plot(), 
       device = "png",
       width = 18, 
       height = 7, 
       units = "in",
       dpi = 300)

### Wild Basil ####
# NA's removed
basil18 <- basil18[complete.cases(basil18), ]

##### 2018 Model Selection ####
# Global model
bas18_mod <- glmmTMB(proportion_het ~ treatment * date, 
                     data = basil18,
                     family = beta_family(),
                     ziformula = ~ treatment * date,
                     contrasts = list(treatment = "contr.sum", 
                                      date = "contr.sum"))

# Model check
summary(bas18_mod)
simulateResiduals(bas18_mod, 
                  plot = T)

# w/o interaction
bas18_mod2 <- glmmTMB(proportion_het ~ treatment + date, 
                      data = basil18,
                      family = beta_family(),
                      ziformula = ~ treatment + date,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))

# Compare
AIC(bas18_mod, 
    bas18_mod2)

# Model check
summary(bas18_mod2)
simulateResiduals(bas18_mod2, 
                  plot = T)

# w/o date
bas18_mod3 <- glmmTMB(proportion_het ~ treatment, 
                      data = basil18,
                      family = beta_family(),
                      ziformula = ~ treatment + date,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))

# Compare
AIC(bas18_mod, 
    bas18_mod2, 
    bas18_mod3) # Best

# Model check
summary(bas18_mod3)
simulateResiduals(bas18_mod3, 
                  plot = T)

###### Best Fit ####
app_summs[[7]] <- bas18_mod3
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment
# Zero inflation:                  ~treatment + date
# Data: basil18
# 
# AIC      BIC   logLik deviance df.resid 
# -8.5      3.8     10.3    -20.5       52 
# 
# 
# Dispersion parameter for beta family (): 5.34 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.1009     0.1328  -8.291  < 2e-16 ***
#   treatment1   -0.3575     0.1257  -2.845  0.00444 ** 
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   -12.291  13140.164  -0.001   0.9993  
# treatment1     11.679  13140.164   0.001   0.9993  
# date1           1.123      0.481   2.334   0.0196 *

###### LRT ####
depo_lrts[[7]] <- lrtest(bas18_mod2, 
                    update(bas18_mod2, 
                           .~. -treatment))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + date
# Model 2: proportion_het ~ date
#   #Df  LogLik Df  Chisq Pr(>Chisq)   
# 1   7 10.3428                        
# 2   6  6.7872 -1 7.1111   0.007661 **

# Determine effect of damage
pairs(emmeans(bas18_mod2, 
              ~ date + treatment, 
              type = "response"))
# contrast         odds.ratio    SE df null t.ratio p.value
# control / damage        0.5 0.128 51    1  -2.703  0.0093
# 
# Results are averaged over the levels of: date 
# Tests are performed on the log odds ratio scale

emmeans(bas18_mod2, 
        ~ treatment + date, 
        type = "response")
# treatment plotpair response      SE  df lower.CL upper.CL
# control   gf1    0.0559 0.00698 225   0.0436   0.0714
# damage    gf1    0.0523 0.00688 225   0.0403   0.0676
# control   gf2    0.0913 0.00976 225   0.0738   0.1124
# damage    gf2    0.0856 0.00925 225   0.0690   0.1056
# control   gf3    0.0764 0.00870 225   0.0609   0.0954
# damage    gf3    0.0715 0.00824 225   0.0569   0.0896
# control   gf4    0.0666 0.00992 225   0.0495   0.0890
# damage    gf4    0.0623 0.00948 225   0.0461   0.0838
# 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 

###### Plot ####
bas18_plot <- plot(ggpredict(bas18_mod3, 
                             ~ date + treatment, 
                             type = "zero_inflated"), 
                   dot.size = 5, 
                   line.size = 2)

# Customize
bas18_plot <- bas18_plot + 
  ggtitle(expression(paste("G) 2018 - ", italic("C. vulgare")))) + 
  labs(x = "Date", 
       y = "Proportion of \n Heterospecific Pollen") + 
  ylim(range(0,.6)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("11 July", 
                              "13 July")) +
  theme(legend.position = "none") +
  annotate(geom = "text", 
           y = .6,
           x = 1.75, 
           label = "Treatment: P = 0.01", 
           size = 8)

# change shape for treatment
bas18_plot <- manual_shape_change(bas18_plot, 17)


#### Whorled Loosestrife ##### 
##### 2019 Model Selection ####
## NOTE: not starting w/ interactions because model won't converge

# Global Model with vetch as resource
loo19_modv <- glmmTMB(proportion_het ~ treatment + date + vetch,
                      data = loose19,
                      family = beta_family(),
                      ziformula = ~ treatment + date,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))

# Model check
summary(loo19_modv) # Not working
simulateResiduals(loo19_modv, 
                  plot = T) # Outlier problem 

# Global Model with milkweed as resource
loo19_modm <- glmmTMB(proportion_het ~ treatment + date + milk_in,
                      data = loose19,
                      family = beta_family(),
                      ziformula = ~ treatment + date,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))

# Model check
summary(loo19_modm)
simulateResiduals(loo19_modm, 
                  plot = T) # Outlier problem

# Compare
AIC(loo19_modv, 
    loo19_modm) # Not sure why these are the same

# Without zero inflation formula
loo19_mod2 <- glmmTMB(proportion_het ~ treatment + date + vetch,
                      data = loose19,
                      family = beta_family(),
                      ziformula = ~ 1,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))



# Model check
summary(loo19_mod2)
simulateResiduals(loo19_mod2, 
                  plot = T) # Outlier problem is gone

# Compare
AIC(loo19_modv, 
    loo19_mod2) # simpler

# Without resources
loo19_mod3 <- glmmTMB(proportion_het ~ treatment + date,
                      data = loose19,
                      family = beta_family(),
                      ziformula = ~ 1,
                      contrasts = list(treatment = "contr.sum", 
                                       date = "contr.sum"))

# Model check
summary(loo19_mod3)
simulateResiduals(loo19_mod3, 
                  plot = T) # Outlier problem is gone

# Compare 
AIC(loo19_mod2, # better right at the border
    loo19_mod3)
#            df       AIC
# loo19_mod2  6 -274.4839
# loo19_mod3  5 -272.5604

###### Best Fit ####
app_summs[[8]] <- loo19_mod2
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment + date + vetch
# Zero inflation:                  ~1
# Data: loose19
# 
# AIC      BIC   logLik deviance df.resid 
# -274.5   -257.9    143.2   -286.5      112 
# 
# 
# Dispersion parameter for beta family (): 18.4 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.49428    0.35734  -4.182 2.89e-05 ***
# treatment1   0.32208    0.13058   2.467   0.0136 *  
# date1       -0.24336    0.11937  -2.039   0.0415 *  
# vetch       -0.02277    0.01146  -1.986   0.0470 *  
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -3.118      0.457  -6.823 8.94e-12 ***

###### LRT #####
# Without resources
depo_lrts[[8]] <- lrtest(loo19_mod2, 
                    update(loo19_mod2, 
                           .~. -vetch))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + date + vetch
# Model 2: proportion_het ~ treatment + date
#   #Df LogLik Df  Chisq Pr(>Chisq)  
# 1   6 143.24                       
# 2   5 141.28 -1 3.9235    0.04762 *

# Without treatment
lrtest(loo19_mod2, 
       update(loo19_mod2,
              .~. -treatment))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + date + vetch
# Model 2: proportion_het ~ date + vetch
#    #Df LogLik Df  Chisq Pr(>Chisq)  
# 1   6 143.24                       
# 2   5 140.40 -1 5.6825    0.01714 *

###### Plot ####
loo19_plot <- plot(ggpredict(loo19_mod2, 
                             ~ vetch + treatment, 
                             type = "fe"), 
                   dot.size = 5, 
                   line.size = 2, 
                   colors = "bw", 
                   show.legend = FALSE)  # this changes the linetype

loo19_plot <- loo19_plot + 
  ggtitle(expression(paste("H) 2019 - ", 
                           italic("L. quadrifolia")))) + 
  labs(x = bquote(atop(paste("No. of ", 
                             italic("Vicia cracca")), 
                       "Flower Racemes")),
       # y = "Proportion of \n Heterospecific Pollen"
       y = ""
  ) + 
  ylim(range(0,.3)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_fill_manual(name = "Treatment", 
                    values = c("control" = cb[4], 
                               "damage" = cb[7]), 
                    labels = c("Control",
                               "Herbivory")) + 
  scale_linetype_manual(name = "Treatment", 
                        values = c("control" = 1, 
                                   "damage" = 2), 
                        labels = c("Control",
                                   "Herbivory")) + 
  theme(legend.position = "none") +
  annotate(geom = "text", 
           y = .27,
           x = 45, 
           label = "Treatment: P = 0.01 \nResources: P = 0.05", 
           size = 8) 

# this gets this in the same dimensions as the other plots
loo19_plot <- as.ggplot(ggplot_gtable(ggplot_build(loo19_plot)))

# get legend for 
loose_legend <-  readRDS("figures/loose_legend.RDS")

##### Two Species Plot ####
two_spec_plots <- (bas18_plot + loo19_plot) /
  loose_legend + 
  theme(plot.title = element_text(size = 30))

# Save
ggsave("figures/clin_lysi_deposition.png", 
       last_plot(), 
       device = "png",
       width = 14, 
       height = 7, 
       units = "in",
       dpi = 300)

#### All Species Plots ####
all_specs_plots <- 
  (vet18_plot + vet19_plot + vet21_plot) / # first row
  (gal18_plot + gal19_plot + gal21_plot) / # second row
  (bas18_plot + loo19_plot + plot_spacer()) / # third row
  (as.ggplot(common_legend) + as.ggplot(loose_legend)) + # legends
  plot_layout(nrow = 4, # specify arrangement
              heights = c(8,8,8,1), 
              widths = c(8,8,8)) & 
  theme(plot.title = element_text(size = 30))

# Save
ggsave("figures/Figure_2.png", 
       last_plot(), 
       device = "png",
       width = 21, 
       height = 25, 
       units = "in",
       dpi = 300)


#### Save Objects for Appendix ####
saveRDS(object = app_summs, 
        file = "data/appendix_summaries.RDS")

saveRDS(object = depo_lrts, 
        file = "data/appendix_depo_lrts.RDS")

# ###### editing legends for ggeffect::plot()
# guides(linetype = guide_legend(override.aes = list(linetype = c(1,1),
#                                                    color = c(cb[4], cb[7]),
#                                                    shape = c(17,19)), 
#                                keywidth = 5))
