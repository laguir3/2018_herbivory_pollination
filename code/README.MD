# Herbivory-Induced Effects on Community-Wide Pollination Services
## Description
This document describes the R scripts with the code used to to perform the statistical analyses in Aguirre and Adler (2023). 

Data collected, processed and analyzed by Luis A. Aguirre.
Updated: 2023-05-10

---

## Pollen deposition 
**File name:**  
- 01_pollen_deposition.R  

This script has three main sections. In each of the sections we analyze how herbivory on a flowering species (*A. syriaca*) affects the deposition of heterospecific pollen in 1) vetch (*V. cracca*), 2) galium (*G. palustre*) and 3) the two species present during only a single field season (*C. vulgare* and *L. quadrifolia*). Furthermore, each analysis section for each species, is broken down by year. We analysed each species data by year because of slight differences in collection protocols during every field season. 

For each year's data, we began with a global model using the `glmmTMB()` function using `family = beta_family()` to account for values bound by 0 and 1. We also included `zi_formula` terms, to account for zero-inflation in the pollen deposition data. We then proceeded to remove terms using model selection. However, we always retained the herbivory treatment (along with any significant terms), as this was our variable of interest. We compared un-reduced and reduced models using Akaike Information Criterion (AIC) and parsimony (i.e., more complex models were deemed better only if their AIC values were 2 units lower than simpler models). We checked the fit of every model using `simulateResiduals()` from the `DHARMa` package. After selecting the best model, we performed *post hoc* likelihood ratio tests for whether models were improved by dropping the `treatment` variable, we used the `lrtest()`. Lastly, we produced plots using functions from the `ggeffects` and `ggplot2` packages. 


## Seed production
**File name:**  
- 02_vetch_seedset.R

This script has only one section. Here, our aim was to determine whether herbivory on a flowering species (*A. syriaca*) affects the production of seeds in a neighboring flower species, vetch (*V. cracca*). Averaged seeds per seed pod were modeled with the `glmmTMB()` function using `family = tweedie()` to account for an abundance 0 values in the data and `dispformula = ~ pods` to account for the increased variance observed as the number of fruits increased. Note that we also attempted to build models with  zero-inflation terms instead of using Tweedie distribution, but zero-inflation models did not perform as well. After declaring the global models, we then proceeded to remove terms using model selection. However, we always retained the herbivory treatment (along with any significant terms), as this was our variable of interest. We compared un-reduced and reduced models using Akaike Information Criterion (AIC) and parsimony (i.e., more complex models were deemed better only if their AIC values were 2 units lower than simpler models). We checked the fit of every model using `simulateResiduals()` from the `DHARMa` package. After selecting the best model, we performed *post hoc* likelihood ratio tests for whether models were improved by dropping the `treatment` variable, we used the `lrtest()`. Lastly, we produced plots using functions from the `ggeffects` and `ggplot2` packages. 

## Pollinia and pollinator composition
**File name:**  
- 03_pollinia.R

This script has four main sections. The first two analysis sections were aimed at determining whether herbivory on the common milkweed (*A. syriaca*) affects how much of this species' pollen was removed by floral visitors. In the first section we performed this analysis at the site level (i.e., pooling all floral visitor together) and in the second section we performed the same analysis but added taxonomic group as a co-factor in the model. For the analyses in first and second sections, we began with global models using the `glmmTMB()` function using `family = nbinom2()` to the account for discrete data that is bound by 0 and is over-dispersed. We also included `zi_formula` terms to account for zero-inflation in the pollinia data.  

The third analysis section (not included in the manuscript) was aimed at determining whether corbiculate bees (i.e., *Bombus* and *Apis*), were more likely to carry pollen in their corbicula in control or herbivory plots. We began this analysis declaring a global model using using `glmmTMB()` function using `family = binomial()` as the response variable was binary (i.e., corbicula pollen present/absent). 

After declaring the global models in the first three sections, we proceeded to remove terms using model selection. However, we always retained the herbivory treatment (along with any significant terms), as this was our variable of interest. We compared un-reduced and reduced models using Akaike Information Criterion (AIC) and parsimony (i.e., more complex models were deemed better only if their AIC values were 2 units lower than simpler models). We checked the fit of every model using `simulateResiduals()` from the `DHARMa` package. After selecting the best model, we performed *post hoc* likelihood ratio tests for whether models were improved by dropping the `treatment` variable, we used the `lrtest()`. Lastly, we produced plots using functions from the `ggeffects` and `ggplot2` packages. 

In the analysis in the fourth section, we examined whether the composition of the observed communities of floral visitors to milkweeds (*A. syriaca*) differed between herbivory and control plots. To do this we conducted Goodness-of-Fit $\chi^2$ tests using the `chisq.test()`. Within this function we specified the observations in the damage plot as the "observed" and the control as the "expected". For this analysis, we chose to pool observations by site and not by date because some taxonomic groups were not observed every day. Lastly, we produced plots summarizing the results using functions from the `ggplot2` package. 

## Appendix 
**File name:**  
- 04_deposition_appendix.R  

This script includes analyses, tables and plots used in the Appendix of Aguirre and Adler (2023). The analysis in the first section, entitled "Effect Spillover", is aimed at determining whether the effect of herbivory found in the previous analysis was the same in the area adjacent to the boundaries of the plot. The first section follows the same procedure explained in the *Pollen deposition* section above.  

In the second section, entitled "Resources", we performed tests to determine whether the amount of resources vary between control and damage plots. Primarily, we rely on the results from the permutation tests using the `symmetry_test()` function from the `coin` package.   

In the next section, entitled "Appendix Tables", we created the majority of the tables used in the appendix, Tables S2-S6. Note that the table contents were derived from the analyses performed in the previous three scripts.   

In the last section, we produced the plots depicting the composition of the observed floral visitor on galium (*G. palustre*) and vetch (*V. cracca*). Note that no analysis accompanies these plots because we did not have enough observations to carry out robust analyses. 
