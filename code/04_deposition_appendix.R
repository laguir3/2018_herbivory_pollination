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
  
  # return
  return(as.ggplot(ggplot_gtable(temp)))
}

# Load datasets w/ inside and outside obs.
# Vicia Cracca
vetch19 <- read.csv("data/2019_vetch_deposition.csv", 
                    header = TRUE)

# Whorled Loosestrife
loose19 <- read.csv("data/2019_loosestrife_deposition.csv",
                    header = TRUE)

# Data Format
# Treatment as factor
vetch19$treatment <- as.factor(vetch19$treatment)
loose19$treatment <- as.factor(loose19$treatment) 

# PlotPair as factor
vetch19$plotpair <- as.factor(vetch19$site)
loose19$plotpair <- as.factor(loose19$site)

### Effect Spillover ####
##### 2019 Vetch #####
###### Damage - In/Out Comparison ####
# Global Model for in/out comparison for damage with vetch resources
spill_dam_modv <- glmmTMB(proportion_het ~ in_out * plotpair + vetch +
                            (1|date),
                          data = filter(vetch19, 
                                        treatment == "damage", # only dates with
                                        date %in% c("2019-07-02", # in and out 
                                                    "2019-07-08", # for damage
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_dam_modv) # no apparent effect
simulateResiduals(spill_dam_modv, 
                  plot = T) # no problems

# Global model with milkweed resources
spill_dam_modm <- glmmTMB(proportion_het ~ in_out * plotpair + milk_in +
                            (1|date),
                          data = filter(vetch19,
                                        treatment == "damage", 
                                        date %in% c("2019-07-02", 
                                                    "2019-07-08",
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_dam_modm)
simulateResiduals(spill_dam_modm, 
                  plot = T) # 

# Compare
AIC(spill_dam_modv, 
    spill_dam_modm) # neither model is good, reduce

# Without random effects
spill_dam_mod2 <- glmmTMB(proportion_het ~ in_out * plotpair + vetch,
                          data = filter(vetch19, 
                                        treatment == "damage", # only dates with
                                        date %in% c("2019-07-02", # in and out 
                                                    "2019-07-08", # for damage
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_dam_mod2)
simulateResiduals(spill_dam_mod2, 
                  plot = T) # 

# Compare
AIC(spill_dam_modv, 
    spill_dam_mod2) # better

# Without interaction
spill_dam_mod3 <- glmmTMB(proportion_het ~ in_out + plotpair + vetch,
                          data = filter(vetch19,
                                        treatment == "damage", 
                                        date %in% c("2019-07-02", 
                                                    "2019-07-08",
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out + plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_dam_mod3)
simulateResiduals(spill_dam_mod3, 
                  plot = T) 

AIC(spill_dam_mod2, 
    spill_dam_mod3) # better

# Without resources
spill_dam_mod4 <- glmmTMB(proportion_het ~ in_out + plotpair,
                          data = filter(vetch19,
                                        treatment == "damage", 
                                        date %in% c("2019-07-02", 
                                                    "2019-07-08",
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out + plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_dam_mod4)
simulateResiduals(spill_dam_mod4, 
                  plot = T) # 

# Compare
AIC(spill_dam_mod3,
    spill_dam_mod4) # better


# Without Plot Pair in main effects
spill_dam_mod5 <- glmmTMB(proportion_het ~ in_out ,
                          data = filter(vetch19,
                                        treatment == "damage", 
                                        date %in% c("2019-07-02", 
                                                    "2019-07-08",
                                                    "2019-07-09")),
                          family = beta_family(),
                          ziformula = ~ in_out,
                          contrasts = list(in_out = "contr.sum"))

# Model check
summary(spill_dam_mod5) # still no discernible effect of in_out
simulateResiduals(spill_dam_mod5, 
                  plot = T) 

# Compare 
AIC(spill_dam_mod4, 
    spill_dam_mod5)
#                df       AIC
# spill_dam_mod4  7 -45.09342
# spill_dam_mod5  5 -44.87531

###### Best Fit ####
# Family: beta  ( logit )
# Formula:          proportion_het ~ in_out
# Zero inflation:                  ~in_out
# Data: 
#   filter(vetch19, treatment == "damage", date %in% c("2019-07-02",  
#                                                      "2019-07-08", 
#                                                      "2019-07-09"))
# 
# AIC      BIC   logLik deviance df.resid 
# -44.9    -34.4     27.4    -54.9       55 
# 
# 
# Dispersion parameter for beta family (): 15.3 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.1633     0.1213 -17.829   <2e-16 ***
# in_out1      -0.1512     0.1083  -1.396    0.163    
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.10059    0.29866  -3.685 0.000229 ***
# in_out1     -0.08899    0.29866  -0.298 0.765725 

###### LRT #### 
lrtest(spill_dam_mod5, 
       update(spill_dam_mod5, 
              .~. -in_out))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ in_out
# Model 2: proportion_het ~ 1
#   #Df LogLik Df  Chisq Pr(>Chisq)
# 1   5 27.438                     
# 2   4 26.481 -1 1.9136     0.1666

# Get estimate for results
ggemmeans(spill_dam_mod4, 
          ~ in_out, 
          type = "zero_inflated")


###### Plot ####
spill_dam_plot <- plot(ggpredict(spill_dam_mod4, 
                                 ~ plotpair + in_out, 
                                 type = "zero_inflated"), 
                       dot.size = 5, 
                       line.size = 2)

# Customize
spill_dam_plot <- spill_dam_plot + 
  ggtitle("Simulated Herbivory Plots") + 
  # xlab("Plot Pair") + 
  xlab("") +
  ylab("Proportion of \n Heterospecific Pollen") +
  ylim(range(0,.15)) +
  theme_classic(base_size = 35) +
  scale_color_manual(name = "Location", 
                     values = c("in" = cb[2], 
                                "out" = cb[6]), 
                     labels = c("Inside Plot",
                                "Outside Plot")) + 
  scale_x_discrete(limits = c("C", 
                              "D")) +
  theme(legend.position = "bottom") +
  annotate(geom = "text", 
           y = .15, 
           x = 1.75, 
           label = "Location: P = 0.16", 
           size = 7)

###### Control - In/Out Comparison ####
# Global Model for in/out comparison  vetch resources
spill_con_modv <- glmmTMB(proportion_het ~ in_out * plotpair + vetch +
                            (1|date),
                          data = filter(vetch19, #only dates DO NOT contain both
                                        treatment == "control", # in and out for
                                        !date %in% c("2019-06-27", # control
                                                     "2019-06-28", 
                                                     "2019-07-05")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_con_modv) # no apparent effect
simulateResiduals(spill_con_modv, 
                  plot = T) # no problems

# Global Model for in/out comparison milk_in resources
spill_con_modm <- glmmTMB(proportion_het ~ in_out * plotpair + milk_in +
                            (1|date),
                          data = filter(vetch19, #only dates DO NOT contain both
                                        treatment == "control", # in and out for
                                        !date %in% c("2019-06-27", # control
                                                     "2019-06-28", 
                                                     "2019-07-05")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_con_modm) # maybe some effect of in and out
simulateResiduals(spill_con_modm, 
                  plot = T) # no problems

# Compare
AIC(spill_con_modv, # essentially same
    spill_con_modm)

# Without random effects
spill_con_mod2 <- glmmTMB(proportion_het ~ in_out * plotpair + vetch,
                          data = filter(vetch19, #only dates DO NOT contain both
                                        treatment == "control", # in and out for
                                        !date %in% c("2019-06-27", # control
                                                     "2019-06-28", 
                                                     "2019-07-05")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_con_mod2) # no apparent effect
simulateResiduals(spill_con_mod2, 
                  plot = T) # no problems

# Compare
AIC(spill_con_modv, 
    spill_con_mod2) # simpler

# Without resources
spill_con_mod3 <- glmmTMB(proportion_het ~ in_out * plotpair,
                          data = filter(vetch19, #only dates DO NOT contain both
                                        treatment == "control", # in and out for
                                        !date %in% c("2019-06-27", # control
                                                     "2019-06-28", 
                                                     "2019-07-05")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_con_mod3) # interaction between in_out and plotpair
simulateResiduals(spill_con_mod3, 
                  plot = T) # no problems

# Compare
AIC(spill_con_mod2, 
    spill_con_mod3) # essentially same, but simpler is better


# Without interaction in fe formula
spill_con_mod4 <- glmmTMB(proportion_het ~ in_out + plotpair,
                          data = filter(vetch19, #only dates DO NOT contain both
                                        treatment == "control", # in and out for
                                        !date %in% c("2019-06-27", # control
                                                     "2019-06-28", 
                                                     "2019-07-05")),
                          family = beta_family(),
                          ziformula = ~ in_out * plotpair,
                          contrasts = list(in_out = "contr.sum",
                                           plotpair = "contr.sum"))

# Model check
summary(spill_con_mod4) # interaction between in_out and plotpair
simulateResiduals(spill_con_mod4, 
                  plot = T) # no problems

# Compare
AIC(spill_con_mod3, # interaction model is better
    spill_con_mod4) 

#               df       AIC
#spill_con_mod3 13 -107.8689
#spill_con_mod4  9 -103.7345

###### Best Fit ####
# Family: beta  ( logit )
# Formula:          proportion_het ~ in_out + plotpair
# Zero inflation:                  ~in_out * plotpair
# Data: filter(vetch19, treatment == "control", !date %in% c("2019-06-27",  
#                                                            "2019-06-28", "2019-07-05"))
# 
# AIC      BIC   logLik deviance df.resid 
# -105.7    -74.8     63.8   -127.7      111 
# 
# 
# Dispersion parameter for beta family (): 21.1 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.55560    0.09906 -25.798  < 2e-16 ***
# in_out1     -0.21391    0.08220  -2.602  0.00926 ** 
# plotpair1        0.32144    0.11211   2.867  0.00414 ** 
# plotpair2       -0.13825    0.11244  -1.230  0.21886    
# 
# Zero-inflation model:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.9340     0.2574  -3.628 0.000286 ***
# in_out1         0.3022     0.2574   1.174 0.240402    
# plotpair1           0.1405     0.3304   0.425 0.670629    
# plotpair2           0.4478     0.2999   1.493 0.135479    
# in_out1:plotpair1   0.2906     0.3304   0.879 0.379134    
# in_out1:plotpair2  -0.6633     0.2999  -2.211 0.027005 *  

###### LRT #####
# Without in_out 
lrtest(spill_con_mod4, 
       update(spill_con_mod4, .~. -in_out))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ in_out + plotpair
# Model 2: proportion_het ~ plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)  
# 1  11 63.846                       
# 2  10 60.536 -1 6.6188    0.01009 *

###### Plot ####
spill_con_plot <- plot(ggpredict(spill_con_mod4, 
                                 ~ plotpair + in_out, 
                                 type = "zero_inflated"), 
                       dot.size = 5, 
                       line.size = 2)

# Customize
spill_con_plot <- spill_con_plot + 
  ggtitle("Control Plots") + 
  xlab("Plot Pair") +
  # ylab("Proportion of \n Heterospecific Pollen") +
  ylab("") +
  ylim(range(0, .15)) +
  theme_classic(base_size = 35) +
  scale_color_manual(name = "Location", 
                     values = c("in" = cb[2], 
                                "out" = cb[6]), 
                     labels = c("Inside Plot",
                                "Outside Plot")) + 
  scale_x_discrete(limits = c("A",
                              "C", 
                              "D")) +
  theme(legend.position = "bottom") +
  annotate(geom = "text", 
           y = 0.15, 
           x = 2.5, 
           label = "Location: P = 0.01", 
           size = 7)


###### In-Only Comparison Between Damage and Control ####
# Global Model with vetch
vet19_in_modv <- glmmTMB(proportion_het ~ treatment * plotpair + vetch +
                           (1|date),
                         data = filter(vetch19,
                                       in_out == "in", 
                                       !date == "2019-07-03"), # no in for damage
                         family = beta_family(),
                         ziformula = ~ treatment * plotpair,
                         contrasts = list(treatment = "contr.sum",
                                          plotpair = "contr.sum"))

# Model check
summary(vet19_in_modv)
simulateResiduals(vet19_in_modv, 
                  plot = T) # A lot of dispersion
testDispersion(vet19_in_modv)
## NOTE: zero inflation estimates, look too suspicious for this model to be good

# Global Model with milk_in
vet19_in_modm <- glmmTMB(proportion_het ~ treatment * plotpair + 
                           log(milk_in) + # Won't converge unless logged
                           (1|date),
                         data = filter(vetch19,
                                       in_out == "in", 
                                       !date == "2019-07-03"), # no in for damage
                         family = beta_family(),
                         ziformula = ~ treatment * plotpair,
                         contrasts = list(treatment = "contr.sum",
                                          plotpair = "contr.sum"))

# Model check
summary(vet19_in_modm)
simulateResiduals(vet19_in_modm, 
                  plot = T) # A lot of dispersion and fail some test of variance
## NOTE: zero inflation estimates, look too suspicious for this model to be good

# Without interaction
vet19_in_mod2 <- glmmTMB(proportion_het ~ treatment + plotpair + vetch +
                           (1|date),
                         data = filter(vetch19,
                                       in_out == "in", 
                                       !date == "2019-07-03"), # no in for damage
                         family = beta_family(),
                         ziformula = ~ treatment + plotpair,
                         contrasts = list(treatment = "contr.sum",
                                          plotpair = "contr.sum"))


# Model check
summary(vet19_in_mod2)
simulateResiduals(vet19_in_mod2, 
                  plot = T) # A lot of dispersion and fail some test of variance
## NOTE: zero nflation estimates, look too suspicious for this model to be good

# Compare 
AIC(vet19_in_modv, 
    vet19_in_mod2) # better

# Without random effects
vet19_in_mod3 <- glmmTMB(proportion_het ~ treatment + plotpair + vetch,
                         data = filter(vetch19,
                                       in_out == "in", 
                                       !date == "2019-07-03"), # no in for damage
                         family = beta_family(),
                         ziformula = ~ treatment + plotpair,
                         contrasts = list(treatment = "contr.sum",
                                          plotpair = "contr.sum"))

# Model check
summary(vet19_in_mod3)
simulateResiduals(vet19_in_mod3, 
                  plot = T) # A lot of dispersion

# Compare
AIC(vet19_in_mod2, 
    vet19_in_mod3) # better
#               df       AIC
# vet19_in_mod2 11 -164.2196
# vet19_in_mod3 10 -166.21961

# Without resources
vet19_in_mod4 <- glmmTMB(proportion_het ~ treatment + plotpair,
                         data = filter(vetch19,
                                       in_out == "in", 
                                       !date == "2019-07-03"), # no in for damage
                         family = beta_family(),
                         ziformula = ~ treatment + plotpair,
                         contrasts = list(treatment = "contr.sum",
                                          plotpair = "contr.sum"))

# Model check
summary(vet19_in_mod4)
simulateResiduals(vet19_in_mod4, 
                  plot = T) # A lot of dispersion

# Compare
AIC(vet19_in_mod3, 
    vet19_in_mod4) # delta is 1.97, but simpler is better
#               df       AIC
# vet19_in_mod3 10 -166.2196
# vet19_in_mod4  9 -168.1827

###### Best Fit #####
# Family: beta  ( logit )
# Formula:          proportion_het ~ treatment + plotpair
# Zero inflation:                  ~treatment + plotpair
# Data: filter(vetch19, in_out == "in", !date == "2019-07-03")
# 
# AIC      BIC   logLik deviance df.resid 
# -168.2   -143.1     93.1   -186.2      111 
# 
# 
# Dispersion parameter for beta family (): 18.2 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.40719    0.09179 -26.225   <2e-16 ***
# treatment1  -0.18760    0.07675  -2.444   0.0145 *  
# plotpair1        0.09911    0.10563   0.938   0.3481    
# plotpair2       -0.21045    0.10489  -2.006   0.0448 *  
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.5185     0.2770  -5.481 4.22e-08 ***
# treatment1    0.4840     0.2410   2.008   0.0446 *  
# plotpair1        -0.5140     0.3874  -1.327   0.1846    
# plotpair2         0.4519     0.3190   1.417   0.1566    

###### LRT #####
lrtest(vet19_in_mod4, 
       update(vet19_in_mod4, .~. -treatment))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ treatment + plotpair
# Model 2: proportion_het ~ plotpair
#   #Df LogLik Df  Chisq Pr(>Chisq)  
# 1   9 93.091                       
# 2   8 90.119 -1 5.9448    0.01476 *

###### Plot ####
vet19_in_plot <- plot(ggpredict(vet19_in_mod4, 
                                ~ plotpair + treatment, 
                                type = "zero_inflated"), 
                      dot.size = 5, 
                      line.size = 2)

# Customize
vet19_in_plot <- vet19_in_plot + 
  ggtitle("") +
  xlab("Plot Pair") +
  # xlab("") +
  ylab("Proportion of \n Heterospecific Pollen") +
  # ylab("") +
  ylim(range(0,.15)) +
  theme_classic(base_size = 20) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("B",
                              "C", 
                              "D")) +
  theme(legend.position = "none") +
  annotate(geom = "text", 
           y = .15, 
           x = 2.5, 
           label = "Treatment: P = 0.02", 
           size = 8)

# Change point shape by treatment 
vet19_in_plot <- manual_shape_change(vet19_in_plot, 17)

vet19_in_plot <- vet19_in_plot / common_legend + 
  plot_layout(heights = c(8,1))

# save
ggsave("figures/vet19_in_plot.png",
       last_plot(),
       device = "png",
       width = 7,
       height = 5,
       units = "in",
       dpi = 300)


# ##### All Effect Spillover Plots #####
# spill_plots <- spill_dam_plot +
#   spill_con_plot  +
#   vet19_in_plot +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 25))
# 
# # Save
# ggsave("figures/spillover_deposition.png",
#        last_plot(),
#        device = "png",
#        width = 18,
#        height = 7,
#        units = "in",
#        dpi = 300)

#### 2019 Whorled Loosestrife ####
# Global_model 
spill_loo_mod <- glmmTMB(proportion_het ~ in_out * vetch + 
                           (1|date), # will not converge with date as fe 
                         data = filter(loose19, # only damage for in/out comp
                                       treatment == "damage"),
                         family = beta_family(),
                         ziformula = ~ 1, # because no inflation in last model
                         contrasts = list(in_out = "contr.sum"))

# Model check
summary(spill_loo_mod)
simulateResiduals(spill_loo_mod, 
                  plot = T) 

# Without random effects
spill_loo_mod2 <- glmmTMB(proportion_het ~ in_out * vetch, 
                          data = filter(loose19, # only damage for in/out comp
                                        treatment == "damage"),
                          family = beta_family(),
                          ziformula = ~ 1, 
                          contrasts = list(in_out = "contr.sum"))

# Model check
summary(spill_loo_mod2)
simulateResiduals(spill_loo_mod2, 
                  plot = T) # outlier problem

# Compare
AIC(spill_loo_mod,
    spill_loo_mod2) # simpler is better 

# Without interaction
spill_loo_mod3 <- glmmTMB(proportion_het ~ in_out + vetch, 
                          data = filter(loose19, # only damage for in/out comp
                                        treatment == "damage"),
                          family = beta_family(),
                          ziformula = ~ 1, 
                          contrasts = list(in_out = "contr.sum"))

# Model check
summary(spill_loo_mod3)
simulateResiduals(spill_loo_mod3, 
                  plot = T) # outlier problem

# Compare
AIC(spill_loo_mod2,
    spill_loo_mod3) # simpler is better 

# Without in_out
spill_loo_mod4 <- glmmTMB(proportion_het ~  vetch, 
                          data = filter(loose19, # only damage for in/out comp
                                        treatment == "damage"),
                          family = beta_family(),
                          ziformula = ~ 1)

# Model check
summary(spill_loo_mod4)
simulateResiduals(spill_loo_mod4, 
                  plot = T) # outlier problem

# Compare
AIC(spill_loo_mod3,
    spill_loo_mod4) # better

# Without resources
spill_loo_mod5 <- glmmTMB(proportion_het ~ 1, 
                          data = filter(loose19, # only damage for in/out comp
                                        treatment == "damage"),
                          family = beta_family(),
                          ziformula = ~ 1, 
                          contrasts = list(in_out = "contr.sum"))

# Model check
summary(spill_loo_mod5)
simulateResiduals(spill_loo_mod5, 
                  plot = T) # outlier problem

# Compare
AIC(spill_loo_mod4,
    spill_loo_mod5) # better

###### Best Fit ####
# Family: beta  ( logit )
# Formula:          proportion_het ~ 1
# Zero inflation:                  ~1
# Data: filter(loose19, treatment == "damage")
# 
# AIC      BIC   logLik deviance df.resid 
# -175.0   -167.8     90.5   -181.0       77 
# 
# 
# Dispersion parameter for beta family (): 14.9 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.24564    0.09771  -22.98   <2e-16 ***
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   -2.944      0.513   -5.74 9.48e-09 ***

###### LRT ####
lrtest(spill_loo_mod3, 
       update(spill_loo_mod3, 
              .~. -in_out))
# Likelihood ratio test
# 
# Model 1: proportion_het ~ in_out + vetch
# Model 2: proportion_het ~ vetch
#   #Df LogLik Df  Chisq Pr(>Chisq)
# 1   5 91.758                     
# 2   4 91.740 -1 0.0373     0.8468

texreg(spill_loo_mod3, 
       file = "figures/spill_loo_mod3.docx")


######## Plot ####
loo19_spill_plot <- plot(ggpredict(spill_loo_mod3, 
                                   ~ vetch + in_out, 
                                   type = "zero_inflated"), 
                         dot.size = 5, 
                         line.size = 2)

# Customize
loo19_spill_plot <- loo19_spill_plot + 
  ggtitle(expression(paste(italic("L. quadrifolia, "), 
                           "Simulated Herbivory Plots"))) + 
  xlab(bquote(atop("Floral Resources (No. of ", 
                   italic("Vicia cracca") ~ "Flower Racemes)"))) +
  ylab("Proportion of \n Heterospecific Pollen") +
  ylim(range(0,.18)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Location", 
                     values = c("in" = cb[2], 
                                "out" = cb[6]), 
                     labels = c("Inside Plot",
                                "Outside Plot")) + 
  theme(legend.position = "bottom") +
  annotate(geom = "text", 
           y = .15, 
           x = 24, 
           label = "Location: P = 0.85 \n Resources: P = 0.11", 
           size = 5.5)

# Save
ggsave("figures/loo19_spill_plot.png",
       last_plot(),
       device = "png",
       width = 10,
       height = 7,
       units = "in",
       dpi = 300)


####### Appendix Table ####
# For loop for making and saving tables
tabs <- vector(mode = "list", 
               length = 8)

# name list objects
names(tabs) <- names(app_summs)

# insert changes to flextable::as_flextable.glmmTMB() to label zero-inflation
# effects: Changes only last till session is closed
trace(what = flextable:::as_flextable.glmmTMB, 
      where = flextable, 
      edit = TRUE)

# In line 10 add: 
# data_t$effect[data_t$component %in% "zi"] <- "Zero-inflation effects"

# In line 26 change to:
# ft <- set_header_labels(ft, term = "Terms", estimate = "Estimate", 
# std.error = "Standard Error", statistic = "z-statistic", p.value = "p-value")

for(i in 1:length(app_summs)){
  # Get flextable
  tabs[[i]] <- flextable::as_flextable(app_summs[[i]])
  
  # Add title at the top
  tabs[[i]] <- add_header_lines(tabs[[i]],
                                values = names(tabs)[i])
  
  # Check Plots
  plot(tabs[[i]])
  
  title = paste0("table_", i)
  
  # Save
  save_as_docx(title = tabs[[i]],
               path = paste0(getwd(), "/tables/table_", i, ".docx"))
}

# Join table for species
# Vetch
read_docx() %>% 
  body_add_flextable(value = tabs[[1]]) %>% 
  body_add_flextable(value = tabs[[2]]) %>% 
  body_add_flextable(value = tabs[[3]]) %>%
  print(target = "tables/appendix_vetch_deposition.docx")

# Galium
read_docx() %>% 
  body_add_flextable(value = tabs[[4]]) %>% 
  body_add_flextable(value = tabs[[5]]) %>% 
  body_add_flextable(value = tabs[[6]]) %>%
  print(target = "tables/appendix_galium_deposition.docx")

read_docx() %>% 
  body_add_flextable(value = tabs[[7]]) %>% 
  body_add_flextable(value = tabs[[8]]) %>% 
  print(target = "tables/appendix_other_deposition.docx")
