library(dbplyr)
library(dplyr)

#map the genes to traits by their snps
#coronary
coronary_input <- read.csv('results/final_results/coronary_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05) %>% 
  select(level, trait, adj_pval)
sig_trait_snp_l0 <- read.csv("results/final_results/coronary_output/level0_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l0$level <- "level0"
gene_snp_l0 <- read.csv("results/final_results/coronary_output/level0_snp_gene.txt", sep="\t") 
level0 <- left_join(sig_trait_snp_l0, gene_snp_l0, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l1 <- read.csv("results/final_results/coronary_output/level1_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l1$level <- "level1"
gene_snp_l1 <- read.csv("results/final_results/coronary_output/level1_snp_gene.txt", sep="\t") 
level1 <- left_join(sig_trait_snp_l1, gene_snp_l1, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l2 <- read.csv("results/final_results/coronary_output/level2_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l2$level <- "level2"
gene_snp_l2 <- read.csv("results/final_results/coronary_output/level2_snp_gene.txt", sep="\t") 
level2 <- left_join(sig_trait_snp_l2, gene_snp_l2, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l3 <- read.csv("results/final_results/coronary_output/level3_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l3$level <- "level3"
gene_snp_l3 <- read.csv("results/final_results/coronary_output/level3_snp_gene.txt", sep="\t") 
level3 <- left_join(sig_trait_snp_l3, gene_snp_l3, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l4 <- read.csv("results/final_results/coronary_output/level4_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l4$level <- "level4"
gene_snp_l4 <- read.csv("results/final_results/coronary_output/level4_snp_gene.txt", sep="\t") 
level4 <- left_join(sig_trait_snp_l4, gene_snp_l4, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
trait_gene_all <- rbind(level0, level1, level2, level3, level4)

coronary_trait_gene <- left_join(coronary_input, trait_gene_all)
write_csv(coronary_trait_gene, "results/final_results/coronary_trait_gene.txt")

#aorta
aorta_input <- read.csv('results/final_results/aorta_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05) %>% 
  select(level, trait, adj_pval)
sig_trait_snp_l0 <- read.csv("results/final_results/aorta_output/level0_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l0$level <- "level0"
gene_snp_l0 <- read.csv("results/final_results/aorta_output/level0_snp_gene.txt", sep="\t") 
level0 <- left_join(sig_trait_snp_l0, gene_snp_l0, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l1 <- read.csv("results/final_results/aorta_output/level1_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l1$level <- "level1"
gene_snp_l1 <- read.csv("results/final_results/aorta_output/level1_snp_gene.txt", sep="\t") 
level1 <- left_join(sig_trait_snp_l1, gene_snp_l1, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l2 <- read.csv("results/final_results/aorta_output/level2_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l2$level <- "level2"
gene_snp_l2 <- read.csv("results/final_results/aorta_output/level2_snp_gene.txt", sep="\t") 
level2 <- left_join(sig_trait_snp_l2, gene_snp_l2, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l3 <- read.csv("results/final_results/aorta_output/level3_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l3$level <- "level3"
gene_snp_l3 <- read.csv("results/final_results/aorta_output/level3_snp_gene.txt", sep="\t") 
level3 <- left_join(sig_trait_snp_l3, gene_snp_l3, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
sig_trait_snp_l4 <- read.csv("results/final_results/aorta_output/level4_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l4$level <- "level4"
gene_snp_l4 <- read.csv("results/final_results/aorta_output/level4_snp_gene.txt", sep="\t") 
level4 <- left_join(sig_trait_snp_l4, gene_snp_l4, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene) %>% 
  unique()
trait_gene_all <- rbind(level0, level1, level2, level3, level4)

aorta_trait_gene <- left_join(aorta_input, trait_gene_all)
write_csv(aorta_trait_gene, "results/final_results/aorta_trait_gene.txt")
