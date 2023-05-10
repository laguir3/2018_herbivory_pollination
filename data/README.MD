---
output:
  html_document: default
  pdf_document: default
---
# DESCRIPTION
Data files associated with Aguirre and Adler (2023). Herbivory-Induced Effects on Community-Wide Pollination Services. In Prep. 

Data collected, processed and analyzed by Luis A. Aguirre.

---
# METADATA

## Deposition data
Files containing the data for analyses performed in R scripts 01_pollen_deposition.R and 04_deposition_appendix.R. 

**List of variables:**  
**pid:** unique plant identifier.  
**treatment:** simulated herbivory treatment, 'damage' or 'control'.  
**date:** date of observation.  
**site:** plot pair identity, named changed to 'plotpair' in scripts.  
**subsample (2018 files only):** a = collected in the morning, b = collected in afternoon.  
**in_out (2019 files only):** collection location, within ('in') or outside ('out') plot borders.   
**conspecific:** number of conspecific pollen grains on stigma.  
**other:** number of pollen from other species on stigma.  
**proportion_cons:** proportion of conspecific pollen grain on stigma (no. conspecific/no. total).  
**proportion_het:** proportion of heterospecific pollen grains (no. heterospecific/no. total).  
**total:** total number of pollen grains on stigma.  
**milk_in (2019 and 2021 files only):** number of *Asclepias syriaca* inflorescences in bloom, within the borders of the study plot.   
**vetch (2019 and 2021 files only):** number of *Vicia cracca* inflorescences in bloom, within the border of the study plot.  
**milk_out (2019 and 2021 files only):** number of *Asclepias syriaca* inflorescences in bloom, within a 5 meter strip outside the plot border.  
**galium (2019 and 2021 files only):** proportion of ground cover with *Galium palustre* plant in bloom.  

**List of files:**   
- 2018_basil_deposition.csv  
- 2018_galium_deposition.csv   
- 2018_vetch_deposition.csv   
- 2019_galium_deposition.csv  
- 2019_loosestrife_deposition.csv  
- 2019_vetch_deposition.csv  
- 2021_galium_deposition.csv  
- 2021_vetch_deposition.csv  

## Vetch seed production data
Files containing the data for analyses performed in R scripts 02_vetch_seedset.R

**List of variables:**  
**vid:** unique plant identifier.   
**treatment:** simulated herbivory treatment, 'damage' or 'control'.   
**site:** plot pair identity, named changed to 'plotpair' in some scripts.    
**pods:** number of seed pods produced.  
**seeds:** averaged seeds per seed pod produced.  
**start (2021 only):**  number of flowers opened when inflorescece was bagged.   

**List of files:**   
- 2018_vetch_seedset.csv  
- 2021_vetch_seedset.csv  

Note: Seedset data was not collected in 2019 because all seed pods where lost due to aphid infestation and mold.

## Pollinia and pollinator data
File containing the data for analyses performed in R scripts 03_pollinia.R  
**List of variables:**  
**bid:** floral visitor unique identifier.  
**site:** plot pair identity, named changed to 'plotpair' in some scripts.  
**treatment:** simulated herbivory treatment, 'damage' or 'control'.   
**date:** date of observation.    
**time:** time of observation.   
**plant_species:** plant species of last observed visit.  
**poll_species:** pollinator (floral visitor) identified to species, if possible.   
**poll_genus:** pollinator (floral visitor) identified to genus, if possible.   
**poll_simple:** simplified taxonomic groups for floral visitors.   
**pollinia:** number of pollinia attached to floral visitor.    
**corbs:** corbicula on floral visitors. '0' if floral visitor did not have corbiculae, '1' if corbicula were found.  

**List of files:**  
- 2018_pollinators.csv  
- 2019_pollinators.csv  
- 2021_pollinators.csv  

## Appendix LRT summaries
Files used in 04_deposition_appendix.R. Files contain summary statistics for the likelihood ratio tests performed for glmm's in scripts 01-03.   

**List of files:**  
- appendix_summaries.RDS  
- appendix_depo_lrts.RDS  
- appendix_seed_lrts.RDS  
- appendix_poll_lrts.RDS  

## Floral resources Data

**List of variables:**  
**treatment:** simulated herbivory treatment, 'damage' or 'control'.   
**date:** date of observation.   
**site:** plot pair identity, named changed to 'plotpair' in some scripts.  
**milk_in (2019 and 2021 files only):** number of *Asclepias syriaca* inflorescences in bloom, within the borders of the study plot.   
**vetch (2019 and 2021 files only):** number of *Vicia cracca* inflorescences in bloom, within the border of the study plot.  
**milk_out (2019 and 2021 files only):** number of *Asclepias syriaca* inflorescences in bloom, within a 5 meter strip outside the plot border.  
**galium (2019 and 2021 files only):** proportion of cover with *Galium palustre* plant in bloom.  

**List of files:**  
- resources.csv  
