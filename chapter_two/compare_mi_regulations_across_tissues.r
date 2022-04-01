library(dbplyr)
library(dplyr)
library(ggplot2)
library(gplots)


#comparing the SNP-gene regulations of myocardial infarction (only enriched in level0) of both tissues
#select "MI", then merge the trait_snps and snp_gene files

#aorta
lv0_trait_snps_aorta <- read.csv('results/final_results/aorta_output/level0_sig_trait_snps.txt', sep="\t") %>%
         filter(trait == "Myocardial infarction")
lv0_snp_gene_aorta <- read.csv('results/final_results/aorta_output/level0_snp_gene.txt', sep="\t")
aorta_mi_snp_gene <- left_join(lv0_trait_snps_aorta, lv0_snp_gene_aorta, by = c('SNPS' = 'snp')) %>% 
  select(SNPS, gene)


#coronary
lv0_trait_snps_coronary <- read.csv('results/final_results/coronary_output/level0_sig_trait_snps.txt', sep="\t") %>%
  filter(trait == "Myocardial infarction")
lv0_snp_gene_coronary <- read.csv('results/final_results/coronary_output/level0_snp_gene.txt', sep="\t")
coronary_mi_snp_gene <- left_join(lv0_trait_snps_coronary, lv0_snp_gene_coronary, by = c('SNPS' = 'snp')) %>% 
  select(SNPS, gene)


intersection <- intersect(aorta_mi_snp_gene, coronary_mi_snp_gene)
unique_aorta <- anti_join(aorta_mi_snp_gene, coronary_mi_snp_gene)
unique_coronary <- anti_join(coronary_mi_snp_gene, aorta_mi_snp_gene)

#venn for eqtl-gene regulation pairs
venn(list(Coronary=1:29, Aorta=4:111))

#venn for eqtl overlaps
venn(list(Coronary=1:27, Aorta=4:78))

#venn for gene overlaps
venn(list(Coronary=1:19, Aorta=3:65))
