library(dbplyr)
library(dplyr)
library(ggpubr)


#coronary bootstraps for enrichment significant, counting the frequency of trait of particular levels out of 10k as the P-values.

og_input <- read.csv("results/coronary/coronary_updated_bs_enrichment/1/significant_enrichment.txt", sep="\t") %>% 
  select(level, trait)
string_dir <- c("results/coronary/coronary_updated_bs_enrichment/")


for (i in 2:10000) {
  add_on_dir <- paste(string_dir, i, "/significant_enrichment.txt", sep="")
  if(!file.exists(add_on_dir)) next
  add_on_input <- read.csv(add_on_dir, sep="\t") %>% 
    select(level, trait)
  og_input <- rbind(og_input, add_on_input)
}


final_freq <- og_input %>% group_by(level, trait) %>% tally() %>% 
  mutate(p_value = n / 10000)

write_csv(final_freq, "results/coronary/coronary_bootstraps_enrichment_pvalue.txt")
coronary_input <- read.csv('results/final_results/coronary_output/significant_enrichment.txt', sep="\t")

#final_freq <- read.csv("results/coronary/coronary_bootstraps_enrichment_pvalue.txt")
coronary_final_enrichment <- left_join(coronary_input, final_freq, keep = TRUE) %>% 
  mutate(level = level.x, trait = trait.x, bootstrap_pval = p_value) %>% 
  select(level, trait, trait_eqtls, adj_pval, bootstrap_pval)
coronary_final_enrichment[is.na(coronary_final_enrichment)] <- 0 

write_csv(coronary_final_enrichment, "results/final_results/coronary_sig_enrichment&bs.txt")


#aorta bootstraps for enrichment significant, counting the frequency of trait of particular levels out of 10k as the P-values.

og_input <- read.csv("results/aorta/aorta_bootstraps_enrichment/1/significant_enrichment.txt", sep="\t") %>% 
  select(level, trait)
string_dir <- c("results/aorta/aorta_bootstraps_enrichment/")


for (i in 2:10000) {
  add_on_dir <- paste(string_dir, i, "/significant_enrichment.txt", sep="")
  if(!file.exists(add_on_dir)) next
  add_on_input <- read.csv(add_on_dir, sep="\t") %>% 
    select(level, trait)
  og_input <- rbind(og_input, add_on_input)
}


final_freq <- og_input %>% group_by(level, trait) %>% tally() %>% 
  mutate(p_value = n / 10000)

write_csv(final_freq, "results/aorta/aorta_bootstraps_enrichment_pvalue.txt")
aorta_input <- read.csv('results/final_results/aorta_output/significant_enrichment.txt', sep="\t")

#final_freq <- read.csv("results/aorta/aorta_bootstraps_enrichment_pvalue.txt")
aorta_final_enrichment <- left_join(aorta_input, final_freq, keep = TRUE) %>% 
  mutate(level = level.x, trait = trait.x, bootstrap_pval = p_value) %>% 
  select(level, trait, trait_eqtls, adj_pval, bootstrap_pval)
aorta_final_enrichment[is.na(aorta_final_enrichment)] <- 0 

write_csv(aorta_final_enrichment, "results/final_results/aorta_sig_enrichment&bs.txt")
