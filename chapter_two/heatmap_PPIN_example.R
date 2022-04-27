library(tidyverse)
library(cowplot)
library(pheatmap)
library(RColorBrewer)



#merge snp + gene + trait, only select gene and trait then use the repeated gene (number of eqtls) as the frequency as gradient"
#map the genes to traits by their snps
#coronary
coronary_input <- read.csv('results/final_results/coronary_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05) %>% 
  select(level, trait)
sig_trait_snp_l0 <- read.csv("results/final_results/coronary_output/level0_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l0$level <- "level0"
gene_snp_l0 <- read.csv("results/final_results/coronary_output/level0_snp_gene.txt", sep="\t") 
level0 <- left_join(sig_trait_snp_l0, gene_snp_l0, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l1 <- read.csv("results/final_results/coronary_output/level1_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l1$level <- "level1"
gene_snp_l1 <- read.csv("results/final_results/coronary_output/level1_snp_gene.txt", sep="\t") 
level1 <- left_join(sig_trait_snp_l1, gene_snp_l1, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l2 <- read.csv("results/final_results/coronary_output/level2_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l2$level <- "level2"
gene_snp_l2 <- read.csv("results/final_results/coronary_output/level2_snp_gene.txt", sep="\t") 
level2 <- left_join(sig_trait_snp_l2, gene_snp_l2, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l3 <- read.csv("results/final_results/coronary_output/level3_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l3$level <- "level3"
gene_snp_l3 <- read.csv("results/final_results/coronary_output/level3_snp_gene.txt", sep="\t") 
level3 <- left_join(sig_trait_snp_l3, gene_snp_l3, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l4 <- read.csv("results/final_results/coronary_output/level4_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l4$level <- "level4"
gene_snp_l4 <- read.csv("results/final_results/coronary_output/level4_snp_gene.txt", sep="\t") 
level4 <- left_join(sig_trait_snp_l4, gene_snp_l4, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
trait_gene_all <- rbind(level0, level1, level2, level3, level4)

coronary_trait_gene_snp <- left_join(coronary_input, trait_gene_all)

#load file
#coronary_trait_gene_input <- read.csv("results/final_results/coronary_trait_gene.txt")

#convert file from rows to columns, add count of genes
coronary_dataset <- table(coronary_trait_gene_snp[c("trait", "gene")]) #this creates the table by grouping these two variables.

#heatmap - level 0
pheatmap (coronary_dataset,
          color = colorRampPalette (c ("white", "red", "dark red","brown", "black")) (31), #keep "white" to make sure zeros are not coloured
          #breaks = c(seq(0,1),seq(2,9),seq(10,19),max(20)),
          cluster_cols = T, # clusters gene columns
          cluster_rows = T,
          cellheight = 9, # specifies size of individual heatmap cells
          cellwidth = 9,
          fontsize = 8,
          filename = "results/final_results/coronary_heatmap(silly).pdf",
          main = "Coronary")

#load aorta file
#aorta_trait_gene_input <- read.csv("results/final_results/aorta_trait_gene.txt")


#aorta
aorta_input <- read.csv('results/final_results/aorta_sig_enrichment&bs.txt', sep=",") %>% 
  filter(bootstrap_pval < 0.05) %>% 
  select(level, trait, adj_pval)
sig_trait_snp_l0 <- read.csv("results/final_results/aorta_output/level0_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l0$level <- "level0"
gene_snp_l0 <- read.csv("results/final_results/aorta_output/level0_snp_gene.txt", sep="\t") 
level0 <- left_join(sig_trait_snp_l0, gene_snp_l0, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l1 <- read.csv("results/final_results/aorta_output/level1_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l1$level <- "level1"
gene_snp_l1 <- read.csv("results/final_results/aorta_output/level1_snp_gene.txt", sep="\t") 
level1 <- left_join(sig_trait_snp_l1, gene_snp_l1, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l2 <- read.csv("results/final_results/aorta_output/level2_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l2$level <- "level2"
gene_snp_l2 <- read.csv("results/final_results/aorta_output/level2_snp_gene.txt", sep="\t") 
level2 <- left_join(sig_trait_snp_l2, gene_snp_l2, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l3 <- read.csv("results/final_results/aorta_output/level3_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l3$level <- "level3"
gene_snp_l3 <- read.csv("results/final_results/aorta_output/level3_snp_gene.txt", sep="\t") 
level3 <- left_join(sig_trait_snp_l3, gene_snp_l3, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
sig_trait_snp_l4 <- read.csv("results/final_results/aorta_output/level4_sig_trait_snps.txt", sep="\t") 
sig_trait_snp_l4$level <- "level4"
gene_snp_l4 <- read.csv("results/final_results/aorta_output/level4_snp_gene.txt", sep="\t") 
level4 <- left_join(sig_trait_snp_l4, gene_snp_l4, by = c("SNPS"="snp")) %>% 
  select(level, trait, gene, SNPS)
trait_gene_all <- rbind(level0, level1, level2, level3, level4)

aorta_trait_gene_snp <- left_join(aorta_input, trait_gene_all)


#convert file from rows to columns, add count of genes
aorta_dataset <- table(aorta_trait_gene_snp[c("trait", "gene")]) #this creates the table by grouping these two variables.


my_palette <- c(colorRampPalette(c("white","red", "dark red","brown", "black","blue","blue","blue","blue","blue","blue","blue"))(n=86))


pheatmap (aorta_dataset,
          color = my_palette, #keep "white" to make sure zeros are not coloured
          #breaks = bk,
          cluster_cols = T, # clusters gene columns
          cluster_rows = T,
          cellheight = 9, # specifies size of individual heatmap cells
          cellwidth = 9,
          main = "Aorta",
          filename = "results/final_results/aorta_heatmap(silly).pdf")



