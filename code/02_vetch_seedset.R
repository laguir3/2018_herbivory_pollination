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
library(AICcmodavg)
library(gfplot)
library(ggeffects) # get estimates
library(cowplot) # join plots
library(patchwork) # join plots


# Set color-blind palette
cb <- c("#000000", # black
        "#E69F00", # orange
        "#56B4E9", # light blue
        "#009E73", # green
        "#F0E442", # yellow
        "#0072B2", # dark blue
        "#D55E00", # red
        "#CC79A7") # pink

# Load common legend for plots 
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


# Load data
# 2018
seed18 <- read.csv("data/2018_vetch_seedset.csv",
                   header = T)

# 2021
seed21 <- read.csv("data/2021_vetch_seedset.csv", 
                   header = T)

# Format - 0's instead of NA for seeds per pod
# 2018 data
for(i in 1:nrow(seed18)){
  if(is.na(seed18$seeds[i])){
    seed18$seeds[i] <- 0
  }
}

# 2021 data 
for(i in 1:nrow(seed21)){
  if(is.na(seed21$seeds[i])){
    seed21$seeds[i] <- 0
  }
}


# Fix row vid 51 in seed21
seed21$treatment[51] <- "control"

# levels
seed18$site <- as.factor(seed18$site)
seed21$site <- as.factor(seed21$site)
seed18$treatment <- as.factor(seed18$treatment)
seed21$treatment <- as.factor(seed21$treatment)


# Create list for Likelihood Ratio Test
seed_lrts <- vector("list", 
                    length = 2)

names(seed_lrts) <- c("2018 Vetch Seed Production", 
                      "2021 Vetch Seed Production")

#### MODELS W/ TWEEDIE GAMMA DISTRIBUTION ####
##### 2018 Model Selection ####
# Global Model
seed18_mod <- glmmTMB(seeds ~ treatment * site, 
                      data = seed18, 
                      family = tweedie(),
                      dispformula =  ~ pods + treatment,
                      contrasts = list(treatment = "contr.sum", 
                                       site = "contr.sum"))

# Diagnostics
summary(seed18_mod)
simulateResiduals(seed18_mod, 
                  plot = T)
bptest(seed18_mod)

# Model w/o interaction
seed18_mod2 <- glmmTMB(seeds ~ treatment + site, 
                       data = seed18, 
                       family = tweedie(),
                       dispformula =  ~ pods + treatment,
                       contrasts = list(treatment = "contr.sum", 
                                        site = "contr.sum"))

# Compare
AIC(seed18_mod, 
    seed18_mod2)


###### Best Fit ####
# Family: tweedie( log )
# Formula:          seeds ~ treatment * site
# Dispersion:             ~pods + treatment
# Data: seed18
# 
# AIC      BIC   logLik deviance df.resid 
# 63.6     76.7    -23.8     47.6       30 
# 
# 
# Conditional model:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       0.843604   0.018925   44.58  < 2e-16 ***
# treatment1        0.027934   0.019084    1.46  0.14326    
# site1            -0.002962   0.018873   -0.16  0.87528    
# treatment1:site1  0.052053   0.018842    2.76  0.00573 ** 
# 
# Dispersion model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   2.4922     0.5197   4.795 1.63e-06 ***
# pods         -1.2819     0.1560  -8.216  < 2e-16 ***
# treatment1   -0.7109     0.2864  -2.482   0.0131 *  

###### LRT #####
seed_lrts[[1]] <- lrtest(seed18_mod, 
                         update(seed18_mod, 
                                .~. -treatment:site))
# Likelihood ratio test
# 
# Model 1: seeds ~ treatment * site
# Model 2: seeds ~ treatment + site
# #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1   8 -23.783                       
# 2   7 -26.595 -1 5.6249    0.01771 *

###### Plot ####
seed18_plot <- plot(ggpredict(seed18_mod, 
                              ~ site + treatment, 
                              type = "fe"),
                    dot.size = 5,  
                    line.size = 2)

seed18_plot <- seed18_plot + 
  ggtitle("2018") + 
  labs(x = "Plot Pair",
       y = "Seeds per Fruit") + 
  ylim(range(0,4)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits=c("A", 
                            "C")) +
  theme(legend.position = "none") + 
  annotate(geom = "text", 
           y = 3.5, 
           x = 1.75, 
           label = "Treatment: P = 0.14 \n Treatment x Plot Pair = P = 0.01", 
           size = 7)

# change shape for treatment
seed18_plot <- manual_shape_change(seed18_plot, 17)

##### 2021 Model Selection ####
seed21_mod <- glmmTMB(seeds ~ treatment + site +
                        treatment:site, 
                      data = seed21, 
                      family = tweedie(),
                      dispformula = ~ pods,
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

# Model check
summary(seed21_mod) # estimates look off, this may be a result of overfitting
simulateResiduals(seed21_mod, 
                  plot = T)


# Remove treatment from dispersion formula
seed21_mod2 <-  glmmTMB(seeds ~ treatment + site, 
                        data = seed21, 
                        family = tweedie(),
                        dispformula = ~ pods,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum"))

# Model check
summary(seed21_mod2)
simulateResiduals(seed21_mod2, 
                  plot = T)
testDispersion(simulateResiduals(seed21_mod2, 
                                 plot = T), 
               alternative = "less")

# Compare
AIC(seed21_mod2, 
    seed21_mod3) 

###### Best Fit ####
# Family: tweedie  ( log )
# Formula:          seeds ~ treatment + site
# Dispersion:             ~pods
# Data: seed21
# 
# AIC      BIC   logLik deviance df.resid 
# 186.6    205.4    -85.3    170.6       69 
# 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.91063    0.10817   8.419   <2e-16 ***
# treatment1   0.08067    0.04965   1.625    0.104    
# site1       -0.07957    0.19410  -0.410    0.682    
# site2       -0.28275    0.27106  -1.043    0.297    
# site3        0.12050    0.11582   1.040    0.298    
# 
# Dispersion model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.80612    0.12356   14.62   <2e-16 ***
# pods        -0.71965    0.06649  -10.82   <2e-16 ***

###### LRT #####
seed_lrts[[2]] <- lrtest(seed21_mod2, 
                         update(seed21_mod2, 
                                .~. -treatment))

###### Plot ####
seed21_plot <- plot(ggpredict(seed21_mod2, 
                              terms = c("site", "treatment"),
                              type = "fe"),
                    dot.size = 5, 
                    line.size = 2)

seed21_plot <- seed21_plot + 
  ggtitle("2021") + 
  labs(x = "Plot Pair",
       y = "Seeds per Fruit") + 
  ylim(range(0,6.5)) +
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
  theme(legend.position = "none")  + 
  annotate(geom = "text", 
           y = 6, 
           x = 3, 
           label = "Treatment: P = 0.10", 
           size = 7)

# change shape for treatment
seed21_plot <- manual_shape_change(seed21_plot, 17)

##### All Plots Together ####
seed_plots <- (seed18_plot + seed21_plot) /
  common_legend + 
  plot_annotation(expression(paste(italic("Vicia cracca "), 
                                   "Seed Production"))) + 
  plot_layout(guides = "collect", 
              height = c(8,1)) &
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 30))

# Save
ggsave("figures/Figure_3.png", 
       last_plot(), 
       device = "png",
       width = 18, 
       height = 7, 
       units = "in",
       dpi = 300)

#### Save Objects for Appendix ####
saveRDS(object = seed_lrts, 
        file = "data/appendix_seed_lrts.RDS")
