#### DESCRIPTION ####
# This is the analysis for the pollinator behavior data in Aguirre et al. 
# Effects of herbivory on community pollination processes. 
# 
# Code by Luis A. Aguirre
# 
#### LOAD ####
# Packages 
library(glmmTMB) 
library(tidyverse) # data wrangling
library(DHARMa) # model diagnostics
library(emmeans)
library(car) 
library(lmtest) # likelihood ratio test
library(AICcmodavg)
library(gfplot) # Tweedie dist. residual plots
library(ggeffects) 
library(effects)
library(patchwork)
library(flextable)
library(EnvStats) # Add sample size to barplots
library(ggplotify)
library(cowplot)

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

# load data 
pollinia18 <- read.csv("data/2018_pollinators.csv", 
                       header = T)
pollinia19 <- read.csv("data/2019_pollinators.csv", 
                       header = T)
pollinia21 <- read.csv("data/2021_pollinators.csv", 
                       header = T)

# keep only polls on milkweed
mpollinia18 <- pollinia18 %>% filter(plant_species == "milkweed")
mpollinia19 <- pollinia19 %>% filter(plant_species == "milkweed")
mpollinia21 <- pollinia21 %>% filter(plant_species == "milkweed")

# Keep only complete cases
mpollinia19 <- mpollinia19[complete.cases(mpollinia19),]
mpollinia19$pollinia <- as.integer(mpollinia19$pollinia)

# Site, treatment, date and poll species as factors
vars <- c("site",
          "treatment", 
          "date", 
          "poll_species", 
          "poll_genus", 
          "poll_simple")

mpollinia18[, vars] <- lapply(mpollinia18[, vars], 
                              as.factor)
mpollinia19[, vars] <- lapply(mpollinia19[, vars], 
                              as.factor)
mpollinia21[, vars] <- lapply(mpollinia21[, vars], 
                              as.factor)

# Create list for likelihood ratio test 
poll_lrts <- vector("list",
                    length = 3)

names(poll_lrts) <- c("2018 Pollinia Loads",
                      "2019 Pollinia Loads",
                      "2021 Pollinia Loads")

#### POLLINIA MODELS W/ NEGATIVE BINOMIAL DISTRIBUTION #####
###### 2018 Model Selection #####
poll18_mod <- glmmTMB(pollinia ~ treatment * site + date,  
                      data = mpollinia18,
                      family = nbinom2(), 
                      ziformula = ~ treatment * site, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum", 
                                       date = "contr.sum"))

poll18_mod_poi <- glmmTMB(pollinia ~ treatment * site + date,  
                          data = mpollinia18, 
                          family = poisson(), 
                          ziformula = ~ treatment * site, 
                          contrasts = list(treatment = "contr.sum",
                                           site = "contr.sum", 
                                           date = "contr.sum"))

# Compare
AIC(poll18_mod, # best
    poll18_mod_poi)

# Model Check
summary(poll18_mod)
simulateResiduals(poll18_mod, 
                  plot = T) # perfect

# drop date
poll18_mod2 <- glmmTMB(pollinia ~ treatment * site,
                       ziformula = ~ treatment * site, 
                       data = mpollinia18, 
                       family = nbinom2(), 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

AIC(poll18_mod, 
    poll18_mod2) #better

# Model check 
summary(poll18_mod2)
simulateResiduals(poll18_mod2, 
                  plot = T)

# drop interaction 
poll18_mod3 <- glmmTMB(pollinia ~ treatment + site,
                       ziformula = ~ treatment * site, 
                       data = mpollinia18, 
                       family = nbinom2(), 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))
# Compare
AIC(poll18_mod2, 
    poll18_mod3) # Best

# Model checks
summary(poll18_mod3)
simulateResiduals(poll18_mod3, 
                  plot = T)

# drop ziformula interacton
poll18_mod4 <- glmmTMB(pollinia ~ treatment + site,
                       ziformula = ~ treatment + site, 
                       data = mpollinia18, 
                       family = nbinom2(), 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll18_mod3, 
    poll18_mod4) # better

# Model Check 
summary(poll18_mod4)
simulateResiduals(poll18_mod4, 
                  plot = T)

# drop ziformula treatment
poll18_mod5 <- glmmTMB(pollinia ~ treatment + site,
                       ziformula = ~ site, 
                       data = mpollinia18, 
                       family = nbinom2(), 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll18_mod4, 
    poll18_mod5) # better

# Model Checks 
summary(poll18_mod5)
simulateResiduals(poll18_mod5, 
                  plot = T)

##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment + site
# Zero inflation:            ~site
# Data: mpollinia18
# 
# AIC      BIC   logLik deviance df.resid 
# 1172.0   1195.3   -580.0   1160.0      352 
# 
# 
# Dispersion parameter for nbinom2 family (): 5.67 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 1.613510   0.059459  27.136  < 2e-16 ***
# treatment1  0.159072   0.053254   2.987  0.00282 ** 
# site1       0.004943   0.059538   0.083  0.93384    
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.4490     0.1201   3.738 0.000186 ***
# site1        -0.2650     0.1193  -2.220 0.026404 *  

##### LRT ####
poll_lrts[[1]] <- lrtest(poll18_mod5, 
                         update(poll18_mod5, 
                                .~. -treatment))
# Likelihood ratio test
# 
# Model 1: pollinia ~ treatment + site
# Model 2: pollinia ~ site
#   #Df LogLik Df  Chisq Pr(>Chisq)   
# 1   6 -580.0                        
# 2   5 -584.3 -1 8.5973   0.003367 **

# for example
ggpredict(poll18_mod5, 
          terms = c("treatment"), 
          type = "zero_inflated")

ggpredict(poll18_mod5, 
          terms = c("treatment"), 
          type = "zero_inflated")

ggeffect(poll18_mod5, 
         terms = "treatment", 
         type = "zero_inflated")

ggemmeans(poll18_mod5, 
          ~ treatment)

##### Plots ####
poll18_plot <- plot(ggpredict(poll18_mod5, 
                              ~ site + treatment, 
                              type = "zero_inflated"), 
                    dot.size = 5, 
                    line.size = 2)

# Customize
poll18_plot <- poll18_plot +
  ggtitle("2018") + 
  labs(# Turn on/off
    # x = "Plot Pair",
    x = "", 
    y = "Number of Pollinia \n per Floral Visitor") + 
  ylim(range(0,5)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits=c("A", "C")) +
  theme(legend.position = "") + 
  annotate(geom = "text", 
           y = 5, 
           x = 1.75, 
           label = "Treatment: P = 0.01", 
           size = 7)

# Change shape by treatment
poll18_plot <- manual_shape_change(poll18_plot, 17)

##### 2019 Model Selection ####
poll19_mod <- glmmTMB(pollinia ~ treatment + site + milk_in + 
                        treatment:site + treatment:milk_in  + (1|date),  
                      data = mpollinia19, 
                      family = nbinom2(),
                      ziformula = ~treatment * site, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

poll19_mod_poi <- glmmTMB(pollinia ~ treatment + site + milk_in + 
                            treatment:site + treatment:milk_in  + (1|date),  
                          data = mpollinia19, 
                          family = poisson(),
                          ziformula = ~ treatment * site, 
                          contrasts = list(treatment = "contr.sum",
                                           site = "contr.sum"))

# Compare 
AIC(poll19_mod, # better
    poll19_mod_poi)

# Model check
summary(poll19_mod)
simulateResiduals(poll19_mod, 
                  plot = T)

# drop date 
poll19_mod2 <- glmmTMB(pollinia ~ treatment + site + milk_in + 
                         treatment:site + treatment:milk_in,  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~treatment * site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod, # better
    poll19_mod2)

# Model check
summary(poll19_mod2)
simulateResiduals(poll19_mod2, 
                  plot = T)

# drop interaction in ziformula
poll19_mod3 <- glmmTMB(pollinia ~ treatment + site + milk_in + 
                         treatment:site + treatment:milk_in  + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~treatment + site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod, 
    poll19_mod3) # better

# Model check
summary(poll19_mod3)
simulateResiduals(poll19_mod3, 
                  plot = T)

# drop treatment:site interaction
poll19_mod4 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         treatment:milk_in  + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ treatment + site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod3, 
    poll19_mod4) # better

# Model check
summary(poll19_mod4)
simulateResiduals(poll19_mod4, 
                  plot = T)

# drop treatment from ziformula
poll19_mod5 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         treatment:milk_in  + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod4, 
    poll19_mod5) # better

# Model check
summary(poll19_mod5)
simulateResiduals(poll19_mod5,
                  plot = T)

# drop site from ziformula
poll19_mod6 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         treatment:milk_in  + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ 1, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod5, # better
    poll19_mod6) 

# Model check
summary(poll19_mod6)
simulateResiduals(poll19_mod6, 
                  plot = T) 

# drop treatment:milk_in interaction
poll19_mod7 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))


# Compare
AIC(poll19_mod5, 
    poll19_mod7) # better

# Model check
summary(poll19_mod7)
simulateResiduals(poll19_mod7, # residual vs predicted a little off
                  plot = T)


# drop site 
poll19_mod8 <- glmmTMB(pollinia ~ treatment + milk_in + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod7, 
    poll19_mod8) # better

# Model check
summary(poll19_mod8)
simulateResiduals(poll19_mod8, # residual vs predicted a little off
                  plot = T)

# drop site from zi formula
poll19_mod9 <- glmmTMB(pollinia ~ treatment + milk_in + (1|date),  
                       data = mpollinia19, 
                       family = nbinom2(),
                       ziformula = ~ 1, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll19_mod8, 
    poll19_mod9) # better

# Model check
summary(poll19_mod9)
simulateResiduals(poll19_mod9, # residual vs predicted a little off
                  plot = T)


##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment + site + milk_in + (1 | date)
# Zero inflation:            ~site
# Data: mpollinia19
# 
# AIC      BIC   logLik deviance df.resid 
# 1995.2   2044.2   -985.6   1971.2      426 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# date   (Intercept) 0.06718  0.2592  
# Number of obs: 438, groups:  date, 7
# 
# Dispersion parameter for nbinom2 family (): 2.54 
# 
# Conditional model:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.5640032  0.1669324   9.369  < 2e-16 ***
# treatment1  -0.0815323  0.0464584  -1.755  0.07927 .  
# site1        0.2210771  0.2342254   0.944  0.34524    
# site2       -0.0588341  0.2047976  -0.287  0.77390    
# site3       -0.0069150  0.1732970  -0.040  0.96817    
# milk_in     -0.0020161  0.0007454  -2.705  0.00683 ** 
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.9228     0.3427  -5.610 2.02e-08 ***
# site1        -1.1560     0.6807  -1.698   0.0895 .  
# site2         0.0332     0.4761   0.070   0.9444    
# site3         0.4837     0.3453   1.401   0.1613    

##### LRT #####
poll_lrts[[2]] <- lrtest(poll19_mod7, 
                         update(poll19_mod7, 
                                .~. -treatment))
# Likelihood ratio test
# 
# Model 1: pollinia ~ treatment + site + milk_in + (1 | date)
# Model 2: pollinia ~ site + milk_in + (1 | date)
#   #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1  12 -985.62                       
# 2  11 -987.16 -1 3.0787    0.07933 .

# For results section
ggpredict(poll19_mod7, 
          terms = c("treatment","site"), 
          type = "zero_inflated")

##### Plots ####
poll19_plot <- plot(ggpredict(poll19_mod7, 
                              terms = c("site","treatment"), 
                              type = "zero_inflated"), 
                    dot.size = 5, 
                    line.size = 2)

# Customize
poll19_plot <- poll19_plot + 
  ggtitle("2019") + 
  labs(# Turn on/off
    x = "Plot Pair",
    # y = "Number of Pollinia \n per Floral Visitor", 
    y = "",
  ) + 
  ylim(range(0,7)) +
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
  theme(legend.position = "") +
  annotate(geom = "text", 
           y = 7, 
           x = 3, 
           label = "Treatment: P = 0.08",
           size = 7)

# Change shape by treatment
poll19_plot <- manual_shape_change(poll19_plot, 17)

# Check effect of resource
plot(ggpredict(poll19_mod9, 
               terms = c("milk_in","treatment"), 
               type = "zero_inflated"))

##### 2021 Model Selection ######
poll21_mod <- glmmTMB(pollinia ~ treatment + site + milk_in +
                        treatment:milk_in  + (1|date),  
                      data = mpollinia21, 
                      family = nbinom2(),
                      ziformula = ~ treatment * site, 
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum"))

poll21_mod_poi <- glmmTMB(pollinia ~ treatment + site + milk_in +
                            treatment:milk_in  + (1|date),  
                          data = mpollinia21, 
                          family = poisson(),
                          ziformula = ~ treatment * site, 
                          contrasts = list(treatment = "contr.sum",
                                           site = "contr.sum"))

# Compare 
AIC(poll21_mod,
    poll21_mod_poi)

# Diagnostic
summary(poll21_mod)
simulateResiduals(poll21_mod, 
                  plot = T)

# drop random effects
poll21_mod2 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         treatment:milk_in,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~ treatment * site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll21_mod, 
    poll21_mod2) # better

# Model Check
summary(poll21_mod2)
simulateResiduals(poll21_mod2, 
                  plot = T)

# Drop interaction in zi
poll21_mod3 <- glmmTMB(pollinia ~ treatment + site + milk_in +
                         treatment:milk_in,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~ treatment + site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll21_mod2, 
    poll21_mod3) # better

# Model check
summary(poll21_mod3)
simulateResiduals(poll21_mod3, 
                  plot = T)

# Drop interaction
poll21_mod4 <- glmmTMB(pollinia ~ treatment + site + milk_in,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~treatment + site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll21_mod3, 
    poll21_mod4) # better

# Model Check
summary(poll21_mod4)
simulateResiduals(poll21_mod4, 
                  plot = T)

# Drop treatment from zi 
poll21_mod5 <- glmmTMB(pollinia ~ treatment + site + milk_in,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~ site, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll21_mod4, 
    poll21_mod5) # better

# Model Check
summary(poll21_mod5)
simulateResiduals(poll21_mod5, 
                  plot = T)

# Drop site from zi 
poll21_mod6 <- glmmTMB(pollinia ~ treatment + site + milk_in,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~ 1, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))

# Compare
AIC(poll21_mod5, 
    poll21_mod6) # better

# Model Check
summary(poll21_mod6)
simulateResiduals(poll21_mod6, 
                  plot = T)

# Drop site or resources
poll21_mod7a <- glmmTMB(pollinia ~ treatment + site,  
                       data = mpollinia21, 
                       family = nbinom2(),
                       ziformula = ~ 1, 
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum"))
poll21_mod7b <- glmmTMB(pollinia ~ treatment + milk_in,  
                        data = mpollinia21, 
                        family = nbinom2(),
                        ziformula = ~ 1, 
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum"))


# Compare
AIC(poll21_mod6, 
    poll21_mod7a, 
    poll21_mod7b) # better

# Model Check
summary(poll21_mod7a)
summary(poll21_mod7b)
simulateResiduals(poll21_mod7a, 
                  plot = T)
simulateResiduals(poll21_mod7b, 
                  plot = T)

##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment + site + milk_in
# Zero inflation:            ~1
# Data: mpollinia21
# 
# AIC      BIC   logLik deviance df.resid 
# 415.2    435.9   -199.6    399.2       90 
# 
# 
# Dispersion parameter for nbinom2 family (): 3.59 
# 
# Conditional model:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.911712   0.573819   3.332 0.000864 ***
# treatment1   0.018568   0.104588   0.178 0.859087    
# site1       -0.286562   0.271103  -1.057 0.290502    
# site2        0.182843   0.393627   0.465 0.642284    
# site3       -0.108856   0.188569  -0.577 0.563755    
# milk_in     -0.006009   0.004853  -1.238 0.215683    
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -1.2959     0.3496  -3.706  0.00021 ***

# Note milk_in is significant when we drop the site factor

##### LRT #####
poll_lrts[[3]] <- lrtest(poll21_mod6, 
                         update(poll21_mod6, 
                                .~. -treatment))
# Likelihood ratio test
# 
# Model 1: pollinia ~ treatment + site + milk_in
# Model 2: pollinia ~ site + milk_in
# #Df  LogLik Df  Chisq Pr(>Chisq)
# 1   8 -199.62                     
# 2   7 -199.64 -1 0.0315     0.8592

# Example for results
ggpredict(poll21_mod6, 
          ~treatment, 
          type = "zero_inflated")

##### Plots ####
poll21_plot <- plot(ggpredict(poll21_mod6, 
                              ~ site + treatment, 
                              type = "zero_inflated"), 
                    dot.size = 5, 
                    line.size = 2)

# Customize
poll21_plot <- poll21_plot + 
  ggtitle("2021") + 
  labs(# Turn on/off
    # x = "Plot Pair",
    x = "",
    # y = "Number of Pollinia \n per Floral Visitor", 
    y = "",
  ) + 
  ylim(range(0,7)) +
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
           y = 7, 
           x = 3, 
           label = "Treatment: P = 0.86",
           size = 7)

poll21_plot <- manual_shape_change(poll21_plot, 17)


##### All Plot Together #####
pollinia_plots <- (poll18_plot + poll19_plot + poll21_plot) /
  common_legend  +
  plot_layout(guides = "collect", 
              height = c(8,1)) & 
  theme(plot.title = element_text(size = 30))

# Save
ggsave("figures/Figure_S1.png",
       last_plot(),
       device = "png",
       width = 18,
       height = 7,
       units = "in",
       dpi = 300)


#### POLLINIA BY GENUS MODELS #####
##### Subset Data for Bombus, Apis and Others ####
# Subset 2018 
keep <- c("Apis", "Bombus", "Other Bees")
bees_poll18 <- mpollinia18 %>% 
  mutate(poll_simple = case_when(poll_simple == keep[1] ~
                                   as.character(keep[1]),
                                 poll_simple == keep[2] ~ 
                                   as.character(keep[2]), 
                                 poll_simple == keep[3] ~ 
                                   as.character(keep[3]),
                                 TRUE ~ "Other Floral Visitors"))

# `poll_simple` as factor
bees_poll18$poll_simple <- as.factor(bees_poll18$poll_simple)

# Subset 2019
bees_poll19 <- mpollinia19 %>% # only apis and bombus due to low No
  mutate(poll_simple = case_when(poll_simple == keep[1] ~
                                   as.character(keep[1]),
                                 poll_simple == keep[2] ~ 
                                   as.character(keep[2]), 
                                 TRUE ~ "Other Floral Visitors"))

# `poll_simple` as factor
bees_poll19$poll_simple <- as.factor(bees_poll19$poll_simple)

# Subset 2021
bees_poll21 <- mpollinia21 %>% # only apis and bombus due to low No
  mutate(poll_simple = case_when(poll_simple == keep[1] ~
                                   as.character(keep[1]),
                                 poll_simple == keep[2] ~ 
                                   as.character(keep[2]), 
                                 TRUE ~ "Other Floral Visitors"))

# `poll_simple` as factor
bees_poll21$poll_simple <- as.factor(bees_poll21$poll_simple)

##### 2018 Model Selection ####
# More simple global model due to convergence issues
apoll18_mod <- glmmTMB(pollinia ~ treatment * poll_simple + date,  
                       data = bees_poll18,
                       family = nbinom2(), 
                       ziformula = ~ treatment + poll_simple,
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum"))

# Model check
simulateResiduals(apoll18_mod, 
                  plot = T)
summary(apoll18_mod)

# Drop treatement from zi
apoll18_mod2 <- glmmTMB(pollinia ~ treatment * poll_simple + date,  
                        data = bees_poll18,
                        family = nbinom2(), 
                        ziformula = ~ poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll18_mod2, 
                  plot = T)
summary(apoll18_mod2)

# Compare
AICc(apoll18_mod)
AICc(apoll18_mod2) # Better

# Drop date
apoll18_mod3 <- glmmTMB(pollinia ~ treatment * poll_simple,  
                        data = bees_poll18,
                        family = nbinom2(), 
                        ziformula = ~ poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll18_mod3, 
                  plot = T)
summary(apoll18_mod3)

# Compare
AICc(apoll18_mod2)
AICc(apoll18_mod3)

# Drop interaction
apoll18_mod4 <- glmmTMB(pollinia ~ treatment + poll_simple,  
                        data = bees_poll18,
                        family = nbinom2(), 
                        ziformula = ~ poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll18_mod4, 
                  plot = T)
summary(apoll18_mod4)

# Compare
AICc(apoll18_mod3)
AICc(apoll18_mod4) # better, done

##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment + poll_simple
# Zero inflation:            ~poll_simple
# Data: bees_poll18
# 
# AIC      BIC   logLik deviance df.resid 
# 1042.1   1080.9   -511.0   1022.1      348 
# 
# 
# Dispersion parameter for nbinom2 family (): 6.78 
# 
# Conditional model:
#                                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       1.69963    0.06047  28.108  < 2e-16 ***
# treatment1                        0.17088    0.05151   3.318 0.000908 ***
# poll_simpleBombus                -0.13635    0.12664  -1.077 0.281646    
# poll_simpleOther Bee             -0.54869    0.24301  -2.258 0.023951 *  
# poll_simpleOther Floral Visitors -0.80266    0.33207  -2.417 0.015642 *  
# 
# Zero-inflation model:
#                                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       -1.6007     0.2664  -6.008 1.87e-09 ***
# poll_simpleBombus                  2.9510     0.3318   8.895  < 2e-16 ***
# poll_simpleOther Bee               2.9633     0.4636   6.392 1.63e-10 ***
# poll_simpleOther Floral Visitors   2.9106     0.5457   5.334 9.60e-08 ***

##### Plots ####
apoll18_plot <- plot(ggpredict(apoll18_mod4, 
                              ~ poll_simple + treatment, 
                              type = "fe.zi"), 
                     dot.size = 5, 
                     line.size = 2)

# Get sample size 
nsizes <- table(bees_poll18$poll_simple)
labs <- c("Apis", 
          "Bombus", 
          "Other Bees", 
          "OFV's")

for(i in 1:length(nsizes)){
  labs[i] <- paste(labs[i],
                   " \n n =",
                   nsizes[i])}

# Customize
apoll18_plot <- apoll18_plot + 
  ggtitle("2018") + 
  labs(# Turn on/off
    # x = "Species",
    x = "",
    y = "Number of Pollinia \n per Floral Visitor",
    # y = "",
  ) + 
  ylim(range(0,7)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) +
  scale_x_discrete(limits = c("Apis", 
                              "Bombus",
                              "Other Bees",
                              "Other Floral Visitors"),
                   labels = c(labs[1:4])) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) + 
  annotate(geom = "text", 
           y = 6.5, 
           x = 3.2, 
           label = "Treatment: \n p-value = 0.01",
           size = 7)

# change point shape by treatment
apoll18_plot <- manual_shape_change(apoll18_plot, 17)

##### 2019 Model Selecition ####
# More simple global model due to convergence issues
apoll19_mod <- glmmTMB(pollinia ~ treatment * poll_simple + site + milk_in +
                         (1|date),  
                       data = bees_poll19,
                       family = nbinom2(), 
                       ziformula = ~ treatment + poll_simple,
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll19_mod, 
                  plot = T)
summary(apoll19_mod)

# drop site
apoll19_mod2 <- glmmTMB(pollinia ~ treatment * poll_simple + milk_in + 
                          (1|date),  
                        data = bees_poll19,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum",
                                         poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll19_mod2, 
                  plot = T)
summary(apoll19_mod2)

# Compare
AICc(apoll19_mod)
AICc(apoll19_mod2)

# Drop random effects
apoll19_mod3 <- glmmTMB(pollinia ~ treatment * poll_simple + milk_in,  
                        data = bees_poll19,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll19_mod3, 
                  plot = T)
summary(apoll19_mod3)

# Compare
AICc(apoll19_mod2) # better, keep random effects
AICc(apoll19_mod3)

# Drop interaction
apoll19_mod4 <- glmmTMB(pollinia ~ treatment + poll_simple + milk_in +
                          (1|date),  
                        data = bees_poll19,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum", 
                                         poll_simple = "contr.sum"))

# Model Check
simulateResiduals(apoll19_mod4, 
                  plot = T)
summary(apoll19_mod4)

# Compare
AICc(apoll19_mod2)
AICc(apoll19_mod4) # better

# Drop poll_simple from zi
apoll19_mod5 <- glmmTMB(pollinia ~ treatment + poll_simple + milk_in +
                          (1|date),  
                        data = bees_poll19,
                        family = nbinom2(), 
                        ziformula = ~ treatment,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll19_mod5, 
                  plot = T)
summary(apoll19_mod5)

# Compare
AICc(apoll19_mod2)
AICc(apoll19_mod4) # best
AICc(apoll19_mod5)

# Drop interaction
apoll19_mod6 <- glmmTMB(pollinia ~ treatment + poll_simple +
                          (1|date),  
                        data = bees_poll19,
                        family = nbinom2(), 
                        ziformula = ~ treatment,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum"))

# Model check
simulateResiduals(apoll19_mod6, 
                  plot = T)
summary(apoll19_mod6)

# Compare
AICc(apoll19_mod4)
AICc(apoll19_mod6) # reducing it makes it worse

##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment + poll_simple + milk_in + (1 | date)
# Zero inflation:            ~treatment + poll_simple
# Data: bees_poll19
# 
# AIC      BIC   logLik deviance df.resid 
# 1854.5   1899.4   -916.3   1832.5      427 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# date   (Intercept) 0.07883  0.2808  
# Number of obs: 438, groups:  date, 7
# 
# Dispersion parameter for nbinom2 family (): 3.08 
# 
# Conditional model:
#                Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.4460769  0.1717897   8.418   <2e-16 ***
# treatment1   -0.0910158  0.0434325  -2.096   0.0361 *  
# poll_simple1  0.1395292  0.0991171   1.408   0.1592    
# poll_simple2  0.2222953  0.1685593   1.319   0.1872    
# milk_in      -0.0016719  0.0006858  -2.438   0.0148 *  
# 
# Zero-inflation model:
#               Estimate Std. Error z value Pr(>|z|)  
# (Intercept)    -6.7045   956.6047  -0.007   0.9944  
# treatment1     -0.4622     0.2177  -2.123   0.0337 *
# poll_simple1  -13.8349  1913.2094  -0.007   0.9942  
# poll_simple2    7.1390   956.6047   0.007   0.9940  

##### Plots ####
apoll19_plot <- plot(ggpredict(apoll19_mod4, 
                               ~ poll_simple + treatment, 
                               type = "fe.zi"), 
                     dot.size = 5, 
                     line.size = 2)

# Get sample size 
nsizes <- table(bees_poll19$poll_simple)
labs <- c("Apis", 
          "Bombus", 
          "OFV's")

for(i in 1:length(nsizes)){
  labs[i] <- paste(labs[i],
                   " \n n =",
                   nsizes[i])}

# Customize
apoll19_plot <- apoll19_plot + 
  ggtitle("2019") + 
  labs(# Turn on/off
    # x = "Species",
    x = "",
    # y = "Number of Pollinia \n per Floral Visitor", 
    y = "",
  ) + 
  ylim(range(0,7)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment", 
                     values = c("control" = cb[4], 
                                "damage" = cb[7]), 
                     labels = c("Control",
                                "Herbivory")) + 
  scale_x_discrete(limits = c("Apis", 
                              "Bombus", 
                              "Other Floral Visitors"),
                   labels = c(labs[1], 
                              labs[2], 
                              labs[3])) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) + 
  annotate(geom = "text", 
           y = 6.5, 
           x = 2.7, 
           label = "Treatment: \n p-value = 0.44",
           size = 7)

# change point shape by treatment
apoll19_plot <- manual_shape_change(apoll19_plot, 17)


##### 2021 Model Selection ####
apoll21_mod <- glmmTMB(pollinia ~ treatment * poll_simple + site + milk_in + 
                         (1|date),  
                       data = bees_poll21,
                       family = nbinom2(), 
                       ziformula = ~ treatment + poll_simple,
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum", 
                                        poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll21_mod, 
                  plot = T)
summary(apoll21_mod) # fail to converge

# Drop site
apoll21_mod2 <- glmmTMB(pollinia ~ treatment * poll_simple + milk_in + 
                          (1|date),  
                        data = bees_poll21,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum", 
                                         poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll21_mod2, 
                  plot = T)
summary(apoll21_mod2) # fail to converge

# Drop random effects
apoll21_mod3 <- glmmTMB(pollinia ~ treatment * poll_simple + milk_in,  
                        data = bees_poll21,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum", 
                                         poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll21_mod3, 
                  plot = T)
summary(apoll21_mod3)

# Compare
AICc(apoll21_mod2)
AICc(apoll21_mod3) 

# Instead drop poll_simple from zi
apoll21_mod4 <- glmmTMB(pollinia ~ treatment * poll_simple + milk_in,  
                        data = bees_poll21,
                        family = nbinom2(), 
                        ziformula = ~ treatment,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum", 
                                         poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll21_mod4, 
                  plot = T)
summary(apoll21_mod4)

# Compare
AICc(apoll21_mod3) # better
AICc(apoll21_mod4) 

# Drop interaction instead 
apoll21_mod5 <- glmmTMB(pollinia ~ treatment + poll_simple + milk_in,  
                        data = bees_poll21,
                        family = nbinom2(), 
                        ziformula = ~ treatment + poll_simple,
                        contrasts = list(treatment = "contr.sum",
                                         site = "contr.sum", 
                                         date = "contr.sum", 
                                         poll_simple = "contr.sum"))

# Model check
simulateResiduals(apoll21_mod5, 
                  plot = T)
summary(apoll21_mod5)

# Compare
AICc(apoll21_mod3) # still better
AICc(apoll21_mod5)

##### Best Fit ####
# Family: nbinom2  ( log )
# Formula:          pollinia ~ treatment * poll_simple + milk_in
# Zero inflation:            ~treatment + poll_simple
# Data: bees_poll21
# 
# AIC      BIC   logLik deviance df.resid 
# 353.8    384.8   -164.9    329.8       86 
# 
# 
# Dispersion parameter for nbinom2 family (): 8.43 
# 
# Conditional model:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              1.794047   0.259649   6.910 4.86e-12 ***
# treatment1              -0.368844   0.179708  -2.052  0.04013 *  
# poll_simple1             0.247958   0.183947   1.348  0.17766    
# poll_simple2             0.163279   0.295719   0.552  0.58085    
# milk_in                 -0.006218   0.001567  -3.969 7.21e-05 ***
# treatment1:poll_simple1  0.445725   0.187308   2.380  0.01733 *  
# treatment1:poll_simple2 -0.789658   0.298908  -2.642  0.00825 ** 
# 
# Zero-inflation model:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)    -6.9709  2996.0463  -0.002    0.998
# treatment1     -0.5537     0.4472  -1.238    0.216
# poll_simple1  -15.3661  5992.0926  -0.003    0.998
# poll_simple2    8.2557  2996.0464   0.003    0.998


##### Plots ####
apoll21_plot <- plot(ggpredict(apoll21_mod3, 
                               ~ poll_simple + treatment, 
                               type = "fe.zi"), 
                     dot.size = 5, 
                     line.size = 2) 

# get sample size 
nsizes <- table(bees_poll21$poll_simple)
labs <- c("Apis", 
          "Bombus", 
          "OFV's")

for(i in 1:length(nsizes)){
  labs[i] <- paste(labs[i],
                   " \n n =",
                   nsizes[i])}

# Customize
apoll21_plot <- apoll21_plot + 
  ggtitle("2021") + 
  labs(# Turn on/off
    # x = "Species",
    x = "",
    # y = "Number of Pollinia \n per Floral Visitor", 
    y = "",
  ) + 
  ylim(range(0, 7)) +
  theme_classic(base_size = 30) +
  scale_color_manual(name = "Treatment",
                     values = c("control" = cb[4],
                                "damage" = cb[7]),
                     labels = c("Control",
                                "Herbivory")) +
  scale_x_discrete(limits = c("Apis", 
                              "Bombus", 
                              "Other Floral Visitors"),
                   labels = c(labs[1], 
                              labs[2], 
                              labs[3])) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) + 
  annotate(geom = "text", 
           y = 6.15, 
           x = 2.5, 
           label = "Treatment: p-value = 0.04 \n Treatment x Taxo. Group: \n p-value = 0.01",
           size = 7)


# change point shape by treatment
apoll21_plot <- manual_shape_change(apoll21_plot, 17)

##### All Plots Together #####
apoll_plots <- (apoll18_plot + apoll19_plot + apoll21_plot) / 
  common_legend +
  plot_layout(guides = "collect", 
              height = c(8,1)) & 
  theme(plot.title = element_text(size = 30))

# Save
ggsave("figures/Figure_S2.png",
       last_plot(),
       device = "png",
       width = 18,
       height = 10,
       units = "in",
       dpi = 300)

### CORBICULA MODELS W/ BINOMIAL DISTRIBUTION ####
##### 2018 Model Selection #####
corb18_mod <- glmmTMB(corbs ~ treatment * poll_simple + site + date,  
                      data = filter(bees_poll18, 
                                    poll_simple == "Apis" | 
                                      poll_simple == "Bombus"),
                      family = binomial(),
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum", 
                                       date = "contr.sum",
                                       poll_simple = "contr.sum"))

# Model Check
summary(corb18_mod)
simulateResiduals(corb18_mod, 
                  plot = T)

# Drop site
corb18_mod2 <- glmmTMB(corbs ~ treatment * poll_simple + date,  
                       data = filter(bees_poll18, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model Check
summary(corb18_mod2)
simulateResiduals(corb18_mod2, 
                  plot = T)

# Compare
AICc(corb18_mod)
AICc(corb18_mod2) # Best

# drop interaction
corb18_mod3 <- glmmTMB(corbs ~ treatment + poll_simple + date,  
                       data = filter(bees_poll18, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model Check
summary(corb18_mod3)
simulateResiduals(corb18_mod3, 
                  plot = T)

# Compare
AICc(corb18_mod2)
AICc(corb18_mod3) # Best

# Drop date
corb18_mod4 <- glmmTMB(corbs ~ treatment + poll_simple,  
                       data = filter(bees_poll18, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model Check
summary(corb18_mod4)
simulateResiduals(corb18_mod4, 
                  plot = T)

# Compare
AICc(corb18_mod3)
AICc(corb18_mod4) # Best


##### 2019 Model Selection ####
corb19_mod <- glmmTMB(corbs ~ treatment * poll_simple + site + scale(milk_in) +
                        (1|date),  
                      data = filter(bees_poll19, 
                                    poll_simple == "Apis" | 
                                      poll_simple == "Bombus"),
                      family = binomial(),
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum", 
                                       date = "contr.sum",
                                       poll_simple = "contr.sum"))

# Model Check
summary(corb19_mod)
simulateResiduals(corb19_mod, 
                  plot = T)

#drop random effect
corb19_mod2 <- glmmTMB(corbs ~ treatment * poll_simple + site + scale(milk_in),  
                      data = filter(bees_poll19, 
                                    poll_simple == "Apis" | 
                                      poll_simple == "Bombus"),
                      family = binomial(),
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum", 
                                       date = "contr.sum",
                                       poll_simple = "contr.sum"))

# Model check
summary(corb19_mod2)
simulateResiduals(corb19_mod2, 
                  plot = T)

# Compare
AICc(corb19_mod)
AICc(corb19_mod2) # Better

# Drop interaction
corb19_mod3 <- glmmTMB(corbs ~ treatment + poll_simple + site + scale(milk_in),  
                       data = filter(bees_poll19, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model Check
summary(corb19_mod3)
simulateResiduals(corb19_mod3, plot = T)

# Compare
AICc(corb19_mod2)
AICc(corb19_mod3) # Better

# Drop site
corb19_mod4 <- glmmTMB(corbs ~ treatment + poll_simple + scale(milk_in),  
                       data = filter(bees_poll19, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model check
summary(corb19_mod4)
simulateResiduals(corb19_mod4, 
                  plot = T)

# Compare
AICc(corb19_mod3)
AICc(corb19_mod4) # Better

# drop milk_in
corb19_mod5 <- glmmTMB(corbs ~ treatment + poll_simple ,  
                       data = filter(bees_poll19, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model check
summary(corb19_mod5)
simulateResiduals(corb19_mod5, 
                  plot = T)

# Compare
AICc(corb19_mod4)
AICc(corb19_mod5) # Better

# Drop milk_in
corb19_mod6 <- glmmTMB(corbs ~ treatment * poll_simple ,  
                       data = filter(bees_poll19, 
                                     poll_simple == "Apis" | 
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum",
                                        poll_simple = "contr.sum"))

# Model check
summary(corb19_mod6)
simulateResiduals(corb19_mod5, 
                  plot = T)

# Compare
AICc(corb19_mod5) # Better
AICc(corb19_mod6) 

##### 2021 Model Selection ####
corb21_mod <- glmmTMB(corbs ~ treatment * poll_simple + milk_in + 
                        (1|date),  
                      data = filter(bees_poll21, 
                                    poll_simple == "Apis" |
                                      poll_simple == "Bombus"),
                      family = binomial(),
                      contrasts = list(treatment = "contr.sum",
                                       site = "contr.sum", 
                                       date = "contr.sum", 
                                       poll_simple = "contr.sum"))

# Model check
simulateResiduals(corb21_mod, 
                  plot = T)
summary(corb21_mod) # fail to converge

# Drop random effect
corb21_mod2 <- glmmTMB(corbs ~ treatment * poll_simple + milk_in,
                       data = filter(bees_poll21, 
                                     poll_simple == "Apis" |
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum", 
                                        poll_simple = "contr.sum"))

# Model check
simulateResiduals(corb21_mod2, 
                  plot = T)
summary(corb21_mod2) # fail to converge

# Compare
AICc(corb21_mod)
AICc(corb21_mod2) # Better

# Drop interaction
corb21_mod3 <- glmmTMB(corbs ~ treatment + poll_simple + milk_in,
                       data = filter(bees_poll21, 
                                     poll_simple == "Apis" |
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum", 
                                        poll_simple = "contr.sum"))

# Model
simulateResiduals(corb21_mod3, 
                  plot = T)
summary(corb21_mod3)

# Compare
AICc(corb21_mod2)
AICc(corb21_mod3) # Better

# Drop milk_in 
corb21_mod4 <- glmmTMB(corbs ~ treatment + poll_simple,
                       data = filter(bees_poll21, 
                                     poll_simple == "Apis" |
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum", 
                                        poll_simple = "contr.sum"))

# Model check
simulateResiduals(corb21_mod4, 
                  plot = T)
summary(corb21_mod4)

# Comparison
AICc(corb21_mod3)
AICc(corb21_mod4) # Better

# Drop taxonomic group 
corb21_mod5 <- glmmTMB(corbs ~ treatment,
                       data = filter(bees_poll21, 
                                     poll_simple == "Apis" |
                                       poll_simple == "Bombus"),
                       family = binomial(),
                       contrasts = list(treatment = "contr.sum",
                                        site = "contr.sum", 
                                        date = "contr.sum", 
                                        poll_simple = "contr.sum"))

# Model check
simulateResiduals(corb21_mod5, 
                  plot = T)
summary(corb21_mod5)

# Compare
AICc(corb21_mod4)
AICc(corb21_mod5) # Better


#### Community Composition ####
ftable(bees_poll18[c("date", 
                     "treatment",
                     "poll_simple")])[[1]]
ftable(bees_poll19[c("treatment", 
                     "poll_simple", 
                     "date")])
ftable(bees_poll21[c("treatment", 
                     "poll_simple", 
                     "date")])

# Rename plot pairs for plotting
levels(bees_poll18$site) <- c("A", "C")
levels(bees_poll19$site) <- c("A", "B", "C", "D")
levels(bees_poll21$site) <- c("A", "B", "C", "D")

##### 2018 Model by Site ####
sites18 <- levels(as.factor(bees_poll18$site))

con_tabs18 <- vector(mode = "list",
                     length = length(sites18))

names(con_tabs18) <- sites18
res18 <- con_tabs18

for(i in 1:length(sites18)){
  # Subset data by date
  temp <- filter(bees_poll18,
                 site == sites18[i])
  
  # Make ftable
  con_tabs18[[i]] <- ftable(temp[c("treatment", 
                                   "poll_simple")])
  
  
  # as.matrix and transpose
  con_tabs18[[i]] <- t(as.matrix(con_tabs18[[i]]))
  
  # Check
  # print(con_tabs18[[i]])
  
  # Fisher exact test
  res18[[i]] <- chisq.test(x = as.data.frame(con_tabs18[[i]])$damage, 
                           p = as.data.frame(con_tabs18[[i]])$control, 
                           rescale.p = T)
  
  # Print 
  print(names(res18)[[i]])
  print(res18[[i]])
}

##### 2019 Model by Site ####
sites19 <- levels(as.factor(bees_poll19$site))

con_tabs19 <- vector(mode = "list",
                     length = length(sites19))

names(con_tabs19) <- sites19
res19 <- con_tabs19

for(i in 1:length(sites19)){
  # Subset data by date
  temp <- filter(bees_poll19,
                 site == sites19[i])
  
  # Make ftable
  con_tabs19[[i]] <- ftable(temp[c("treatment", 
                                   "poll_simple")])
  
  
  # as.matrix and transpose
  con_tabs19[[i]] <- t(as.matrix(con_tabs19[[i]]))
  
  # Check
  # print(con_tabs19[[i]])
  
  # Fisher exact test
  res19[[i]] <- chisq.test(x = as.data.frame(con_tabs19[[i]])$damage, 
                           p = as.data.frame(con_tabs19[[i]])$control, 
                           rescale.p = T)
  
  # Print Results
  print(res19[[i]])
}

##### 2021 Model by Site ####
sites21 <- levels(as.factor(bees_poll21$site))

con_tabs21 <- vector(mode = "list",
                     length = length(sites21))

names(con_tabs21) <- sites21
res21 <- con_tabs21

for(i in 1:length(sites21)){
  # Subset data by date
  temp <- filter(bees_poll21,
                 site == sites21[i])
  
  # Make ftable
  con_tabs21[[i]] <- ftable(temp[c("treatment", 
                                   "poll_simple")])
  
  
  # as.matrix and transpose
  con_tabs21[[i]] <- t(as.matrix(con_tabs21[[i]]))
  
  # Check
  # print(con_tabs21[[i]])
  
  # Fisher exact test
  res21[[i]] <- chisq.test(x = as.data.frame(con_tabs21[[i]])$damage, 
                           p = as.data.frame(con_tabs21[[i]])$control,
                           rescale.p = T)
  
  # Print Results
  print(sites21[i])
  print(res21[[i]])
}

##### Create Flextable for Supplemental Part ####
# create data frame for values
chi_tab <- as.data.frame(matrix(nrow = 10, 
                             ncol = 5))
names(chi_tab) <- c("Year", 
                    "Plot Pair", 
                    "chi", 
                    "Degrees of Freedom", 
                    "P")

# join lists
res_all <- c(res18, res19, res21)
names(res_all) <- c("A", "C", 
                    "A", "B", "C", "D", 
                    "A", "B", "C", "D")

# insert year
chi_tab$Year[1:2] <- "2018"
chi_tab$Year[3:6] <- "2019"
chi_tab$Year[7:10] <- "2021"

# for loop to get remaining values 
for(i in 1:length(res_all)){
  
  # get plot pair
  chi_tab$`Plot Pair`[i] <- names(res_all)[[i]]
  
  # get chi-statistic
  chi_tab$chi[i] <- res_all[[i]]$statistic
  
  # get degrees of freedom
  chi_tab$`Degrees of Freedom`[i] <- res_all[[i]]$parameter
  
  # get p-values 
  chi_tab$P[i] <- res_all[[i]]$p.value
}

# round to 3 decimal points
chi_tab$chi <- round(chi_tab$chi, 
                     digits = 3)

chi_tab$P <- round(chi_tab$P, 
                   digits = 3)

# create flextable
chi_ft <- flextable(chi_tab)

# bold column names
chi_ft <- bold(chi_ft,
               part = "header")

# bold statistically significant plot pairs
chi_ft <- bold(chi_ft, 
               ~ P < 0.05)

# italicize statistically significant plot
chi_ft <- italic(chi_ft, 
                 ~ P < 0.10 & P > 0.05)

# change P to P-value
chi_ft <- compose(chi_ft,
                  part = "header",
                  j = "P", 
                  i = 1,
                  as_paragraph("P-value"))

# change chi to X-statistic
chi_ft <- compose(chi_ft, 
                  part = "header", 
                  j = "chi", 
                  i = 1, 
                  as_paragraph(paste("\u03C7", 
                                     "-statistic")))

# Widen x-stat column to fit
chi_ft <- width(chi_ft, j = 3, width = .9)

# final format
chi_ft <- theme_booktabs(chi_ft)

save_as_docx(title = chi_ft,
             path = paste0(getwd(), 
                           "/tables/table_chi_statistic.docx"))


##### Goodness of fit by Year #####
# 2018
con_tab18 <- ftable(bees_poll18[c("treatment", "poll_simple")])
con_tab18 <- t(as.matrix(con_tab18))
chi_18 <- chisq.test(x = as.data.frame(con_tab18)$damage, 
                     p = as.data.frame(con_tab18)$control,
                     rescale.p = T)

con_tab19 <- ftable(bees_poll19[c("treatment", "poll_simple")])
con_tab19 <- t(as.matrix(con_tab19))
chi_19 <- chisq.test(x = as.data.frame(con_tab19)$damage, 
                     p = as.data.frame(con_tab19)$control, 
                     rescale.p = T)

con_tab21 <- ftable(bees_poll21[c("treatment", "poll_simple")])
con_tab21 <- t(as.matrix(con_tab21))
chi_21 <- chisq.test(x = as.data.frame(con_tab21)$damage, 
                     p = as.data.frame(con_tab21)$control, 
                     rescale.p = T)

#### Plots #####
# Create variables for easy plots #
bees_poll18 <- unite(bees_poll18, 
                     treatxsite, 
                     site, treatment,
                     sep = "-", 
                     remove = FALSE)

bees_poll19 <- unite(bees_poll19, 
                     treatxsite, 
                     site,
                     treatment, 
                     sep = "-", 
                     remove = FALSE)

bees_poll21 <- unite(bees_poll21, 
                     treatxsite,
                     site, 
                     treatment, 
                     sep = "-", 
                     remove = FALSE)



## 2018
# Create text df for significant
ann_text1 <- tribble(
  ~site, ~treatment, ~poll_simple,
  #--|--|----
  "C", 1.5, 1.05 
) 

# get levels
levs <- levels(bees_poll18$poll_simple)

# Plot
comp18 <- ggplot(data = bees_poll18, 
                 aes(x = treatment)) + 
  geom_bar(aes(fill = poll_simple),
           position = "fill") + 
  facet_grid(.~site) + 
  ggtitle("2018") + 
  labs(#turn on/off
    # x = "Treatment"
    x = "", 
    y = "Proportion"
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(0, 1.15)) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) +
  scale_fill_manual(name = "Floral Visitor Group", 
                    values = c("Apis" = cb[6], 
                               "Bombus" = cb[2], 
                               "Other Bees" = cb[1], 
                               "Other Floral Visitors" = cb[8])) + 
  scale_x_discrete(labels = c("CON", 
                              "HERB")) +
  geom_text(data = ann_text1, 
            aes(y = poll_simple), 
            label = "*", 
            size = 20)

## 2019
ann_text1 <- tribble(
  ~site, ~treatment, ~poll_simple,
  #--|--|----
  "A", 1.5, 1.05 
) 

ann_text2 <- ann_text1
ann_text2[1] <- "B"

# same levels as 2018 plot for common legend
levels(bees_poll19$poll_simple) <- c(levels(bees_poll19$poll_simple), "Other Bees")
bees_poll19$poll_simple <- factor(bees_poll19$poll_simple, 
                                  levels = c("Apis", "Bombus", "Other Bees", 
                                             "Other Floral Visitors"))

# Plot
comp19 <- ggplot(data = bees_poll19, 
                 aes(x = treatment)) + 
  geom_bar(aes(fill = poll_simple), 
           position = "fill") + 
  facet_grid(.~site) + 
  ggtitle("2019") + 
  labs(#turn on/off
    x = "Treatment",
    # x = "", 
    # y = "Proportion", 
    y = ""
  ) +
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(0, 1.15)) +
  theme_classic(base_size = 30) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) +
  scale_fill_manual(name = "Floral Visitor Group", 
                    values = c("Apis" = cb[6], 
                               "Bombus" = cb[2], 
                               "Other Bees" = cb[1], 
                               "Other Floral Visitors" = cb[8]),
                    drop = FALSE) + 
  scale_x_discrete(labels = c("CON", 
                              "HERB")) +
  geom_text(data = ann_text1, 
            aes(y = poll_simple), 
            label = expression("\u2020"), 
            size = 12) +
  geom_text(data = ann_text2, 
            aes(y = poll_simple), 
            label = "*", 
            size = 20)


## 2021
# Create text df for significant
ann_text1 <- tribble(
  ~site, ~treatment, ~poll_simple,
  #--|--|----
  "A", 1.5, 1.05 
) 

ann_text2 <- ann_text1
ann_text2[1] <- "B"

# same levels as 2018 plot for common legend
levels(bees_poll21$poll_simple) <- c(levels(bees_poll21$poll_simple), "Other Bees")
bees_poll21$poll_simple <- factor(bees_poll21$poll_simple, 
                                  levels = c("Apis", "Bombus", "Other Bees", 
                                             "Other Floral Visitors"))


# Plot
comp21 <- ggplot(data = bees_poll21, 
                 aes(x = treatment)) + 
  geom_bar(aes(fill = poll_simple), 
           position = "fill") +
  facet_grid(.~site)  + 
  ggtitle("2021") + 
  labs(#turn on/off
    # x = "Treatment"
    x = "", 
    # y = "Proportion", 
    y = "" 
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(0, 1.15)) +
  theme_classic(base_size = 30) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust = 1)) +
  scale_fill_manual(name = "Floral Visitor Group", 
                    values = c("Apis" = cb[6], 
                               "Bombus" = cb[2], 
                               "Other Bees" = cb[1], 
                               "Other Floral Visitors" = cb[8]), 
                    drop = FALSE) + 
  scale_x_discrete(labels = c("CON", 
                              "HERB")) + 
  geom_text(data = ann_text1, 
            aes(y = poll_simple), 
            label = "*", 
            size = 20) +
  geom_text(data = ann_text2, 
            aes(y = poll_simple), 
            label = "*", 
            size = 20)

##### All Plots together ####
(comp_plots <- comp18 + comp19 + comp21 +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom", 
         plot.title = element_text(size = 30))
)

# Save
ggsave("figures/Figure_1.png",
       last_plot(),
       device = "png",
       width = 18,
       height = 7,
       units = "in",
       dpi = 300)

##### Vetch Floral Visitors #####
# keep only polls on vetch
vpollinia18 <- pollinia18 %>% 
  filter(plant_species == "vetch")
vpollinia19 <- pollinia19 %>% 
  filter(plant_species == "vetch")
vpollinia21 <- pollinia21 %>% 
  filter(plant_species == "vetch")

# Format Taxonomic Group Labels
# '18 Make broader "Other Insect" group
vpollinia18$poll_genus <- as.factor(vpollinia18$poll_genus)
levels(vpollinia18$poll_genus)[4:6] <- "OFV's"

# '19 Syrphus to "Other Insect" group
vpollinia19$poll_genus <- as.factor(vpollinia19$poll_genus)
levels(vpollinia19$poll_genus)[5:6] <- "OFV's" 

# '21 Syrphidae to "Other Insect"
vpollinia21$poll_genus <- as.factor(vpollinia21$poll_genus)
levels(vpollinia21$poll_genus)[c(6:7)] <- "OFV's"

vpollinia <- cbind(vpollinia18, vpollinia19, vpollinia19)


vetch_leg <- c("Andrena" = cb[1], 
               "Bombus" = cb[2], 
               "OFV's" = cb[3],
               "Agapostemon" = cb[4],
               "Apis" = cb[5],
               "Lepidoptera" = cb[6],
               "Megachile" = cb[7],
               "Osmia" = cb[8])

# Plot proportions
# 2018
vcomp18 <- ggplot(data = vpollinia18, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_genus), 
           position = "fill") + 
  ggtitle("2018") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = ""
    tag = "A"
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30) +
  scale_x_discrete(labels = c("A", "C")) +
  scale_fill_manual(name = "Taxonomic Group", # Include all for single legend
                    values = vetch_leg,
                    limits = names(vetch_leg)) + # includes nonpresent vars.
  annotate(geom = "text", 
           label = paste("n =", 
                         as.character(table(vpollinia18$site))), 
           y = -0.05, 
           x = c(1:2), 
           size = 6)


# 2019
vcomp19 <- ggplot(data = vpollinia19, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_genus), 
           position = "fill") +
  ggtitle("2019") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = "" 
    tag = "B"
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30) +
  scale_x_discrete(labels = "C") +
  scale_fill_manual(name = "Taxonomic Group", # Include all for single legend
                    values = vetch_leg,
                    limits = names(vetch_leg)) +
  annotate(geom = "text",
           label = paste("n =",
                         as.character(table(vpollinia19$site))),
           y = -0.05,
           x = 1,
           size = 6) 


# 2021
vcomp21 <- ggplot(data = vpollinia21, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_genus), 
           position = "fill")   + 
  ggtitle("2021") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = "" 
    tag = "C"
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30)  +
  scale_x_discrete(labels = c("A", 
                              "B", 
                              "C", 
                              "D")) +
  scale_fill_manual(name = "Taxonomic Group", # Include all for single legend
                    values = vetch_leg,
                    limits = names(vetch_leg)) +
  annotate(geom = "text", 
           label = paste("n =", 
                         as.character(table(vpollinia21$site))), 
           y = -0.05, 
           x = c(1:4), 
           size = 6) 

# Collect plots
vetch_comps <- vcomp18 +
  vcomp19 + 
  vcomp21 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 30))


# Save
ggsave("figures/Figure_S3.png",
       last_plot(),
       device = "png",
       width = 18,
       height = 7,
       units = "in",
       dpi = 300)

##### Galium Floral Visitors #####
# keep only polls on galium
gpollinia18 <- pollinia18 %>% 
  filter(plant_species == "galium")
gpollinia19 <- pollinia19 %>% 
  filter(plant_species == "galium")
gpollinia21 <- pollinia21 %>% 
  filter(plant_species == "galium")

# Change taxonomic group labels
for(i in 1:nrow(gpollinia18)){
  if(gpollinia18$poll_simple[i] == "Bombus" ||
     gpollinia18$poll_simple[i] == "Wasp"){
    gpollinia18$poll_simple[i] <- as.character("Hymenoptera")
  }
  if(gpollinia18$poll_simple[i] == "Syrphus"){
    gpollinia18$poll_simple[i] <- "Diptera"
  }
  if(gpollinia18$poll_simple[i] == "Other"){
    gpollinia18$poll_simple[i] <- "OFV's"
  }
}

for(i in 1:nrow(gpollinia19)){
  if(gpollinia19$poll_simple[i] == "Lepidoptera" | 
     gpollinia19$poll_simple[i] == "Other Insect"){
    gpollinia19$poll_simple[i] <- "OFV's"
  }
  if(!(gpollinia19$poll_simple[i] == "Diptera" ||
       gpollinia19$poll_simple[i] == "OFV's")){
    gpollinia19$poll_simple[i] <- as.character("Hymenoptera")
  }
}

for(i in 1:nrow(gpollinia21)){
  if(gpollinia21$poll_simple[i] == "Lepidoptera" |
     gpollinia21$poll_simple[i] == "Other Insect"){
    gpollinia21$poll_simple[i] <- "OFV's"
  }
  if(!(gpollinia21$poll_simple[i] == "Diptera" ||
       gpollinia21$poll_simple[i] == "OFV's")){
    gpollinia21$poll_simple[i] <- as.character("Hymenoptera")
  }
}

# Plot proportions
# 2018
gcomp18 <- ggplot(data = gpollinia18, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_simple), 
           position = "fill") +
  ggtitle("2018 ") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = "",
    tag = "A"
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30) +
  scale_fill_manual(name = "Taxonomic Group", 
                    values = c("Diptera" = cb[1], 
                               "Hymenoptera" = cb[2], 
                               "OFV's" = cb[3])) +
  scale_x_discrete(labels = c("A", 
                              "C")) +
  annotate(geom = "text", 
           label = paste("n =", 
                         as.character(table(gpollinia18$site))), 
           y = -0.05, 
           x = c(1:2), 
           size = 6)


# 2019
gcomp19 <- ggplot(data = gpollinia19, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_simple), 
           position = "fill") +
  ggtitle("2019") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = "",
    tag = "B" 
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30) +
  scale_x_discrete(labels = c("A", 
                              "B", 
                              "C", 
                              "D")) +
  scale_fill_manual(name = "Taxonomic Group", 
                    values = c("Diptera" = cb[1], 
                               "Hymenoptera" = cb[2], 
                               "OFV's" = cb[3])) +
  annotate(geom = "text", 
           label = paste("n =", 
                         as.character(table(gpollinia19$site))), 
           y = -0.05, 
           x = c(1:4), 
           size = 6)

# 2021
gcomp21 <- ggplot(data = gpollinia21, 
                  aes(x = site)) + 
  geom_bar(aes(fill = poll_simple), 
           position = "fill")   + 
  ggtitle("2021") + 
  labs(#turn on/off
    x = "Plot Pair",
    # x = "", 
    y = "Proportion",
    # y = "",
    tag = "C" 
  ) + 
  scale_y_continuous(breaks = seq(0.0, 1, 0.25),
                     limits = c(-0.05, 1)) +
  theme_classic(base_size = 30) +
  scale_x_discrete(labels = c("A", 
                              "B", 
                              "C", 
                              "D")) +
  scale_fill_manual(name = "Taxonomic Group", 
                    values = c("Diptera" = cb[1], 
                               "Hymenoptera" = cb[2], 
                               "OFV's" = cb[3])) +
  annotate(geom = "text", 
           label = paste("n =", 
                         as.character(table(gpollinia21$site))), 
           y = -0.05, 
           x = c(1:4), 
           size = 6)


galium_comps <- gcomp18 + 
  gcomp19 + 
  gcomp21 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 30))

# Save
ggsave("figures/Figure_S4.png",
       last_plot(),
       device = "png",
       width = 18,
       height = 7,
       units = "in",
       dpi = 300)


#### Save Objects for Appendix Tables ####
saveRDS(object = poll_lrts,
        file = "data/appendix_poll_lrts.RDS")
