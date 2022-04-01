library(tidyverse)
library(ggplot2)
library(wesanderson)
library(dbplyr)
library(dplyr)
library(ggpubr)

#correlation plot

co_eq_input <- as_tibble(read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE)) %>% 
  dplyr::select(gencode_id, expression, interaction_type, tissue)
co_eq_input <- distinct(co_eq_input) %>% mutate(gene_id=gsub("\\..*","",gencode_id))


co_eq_input <- co_eq_input %>% 
  mutate(
    log_exp = log(expression),
    log_exp = case_when(
      is.infinite(log_exp) ~ 0,
      TRUE ~ as.numeric(log_exp)
    )
  )

loeuf_input <- read.csv('gnomad.v2.1.1.lof_metrics.by_gene.txt',sep='\t',header=TRUE)
co_corre_input <- merge(co_eq_input, loeuf_input, by.x="gene_id", by.y="gene_id", how='INNER') %>% 
  dplyr::select(tissue, interaction_type, log_exp, oe_lof_upper) %>% 
  filter(tissue == 'Artery_Aorta' | tissue == 'Artery_Coronary' | tissue == 'Artery_Tibial')

ggplot(co_corre_input, aes(oe_lof_upper, log_exp, shape=tissue, colour=tissue, fill=tissue)) +
  geom_smooth(method=lm) +
  geom_point(size=1) +
  theme_classic() +
  facet_wrap(vars(interaction_type))+
  theme(legend.position = 'bottom', text = element_text(size = 14)) +
  labs(x="LOEUF Score", y="(Log) Median Expression (TPM)")+ 
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1")) + 
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling1"))+ 
  stat_cor(method='pearson')

loeuf_input %>% filter(gene== 'NME1' | gene== 'NME1-NME2' | gene== 'NME6' | gene== 'NME7') %>% select(gene, oe_lof_upper)

