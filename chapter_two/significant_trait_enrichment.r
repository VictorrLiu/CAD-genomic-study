library(tidyverse)
library(ggplot2)
library(dbplyr)
library(dplyr)
library(ggpubr)
library(readxl)
library(forcats)

#present in scattered bubble plots
#coronary, all (levels) of the significant enrichment traits out put from the (tissue gene map) CoDeS3D algorithm

coronary_input <- read.csv('results/final_results/coronary_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05)

coronary_input <- coronary_input %>% 
  mutate(
    n_log_p = -log(adj_pval)
  ) %>% 
  mutate(trait = fct_reorder(trait, n_log_p)) 



ggplot(coronary_input, aes(n_log_p, trait, colour=factor(level), size=trait_eqtls)) + 
  theme_minimal() +
  geom_point(alpha=0.5) +
  scale_colour_viridis_d() +
  labs(x="-Log (Adj. p-value)", y="Siginificanly Enriched Traits")  


#aorta

aorta_input <- read.csv('results/final_results/aorta_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05)

aorta_input <- aorta_input %>% 
  mutate(
    n_log_p = -log(adj_pval)
  ) %>% 
  mutate(trait = fct_reorder(trait, n_log_p)) 



ggplot(aorta_input, aes(n_log_p, trait, colour=factor(level), size=trait_eqtls)) + 
  theme_minimal() +
  geom_point(alpha=0.5) +
  scale_colour_viridis_d() +
  labs(x="-Log (Adj. p-value)", y="Siginificanly Enriched Traits")  