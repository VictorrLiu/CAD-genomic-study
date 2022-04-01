library(tidyverse)
library(ggplot2)
library(wesanderson)
library(dbplyr)
library(dplyr)
library(ggpubr)

#violin plot with the expression data from the CoDeS3d output, grouped by interaction types and colored with the Wes Anderson colors palette 

co_eq_input <- as_tibble(read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE))
co_eq_input <- co_eq_input %>% dplyr::select(gene, expression, interaction_type, tissue)

co_eq_input <- distinct(co_eq_input) 
co_eq_input$present=1

co_eq_input <- co_eq_input %>% pivot_wider(values_from = present, names_from = interaction_type) %>%
  mutate(interaction = case_when(
    Cis == 1 & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans',
    Cis == 1 & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Cis & Trans-Interchromosomal', 
    Cis == 1 & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans-Intrachromosomal',
    is.na(Cis) & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Trans',
    Cis == 1 & is.na(`Trans-interchromosomal`) & is.na(`Trans-intrachromosomal`) ~ 'Cis',
    is.na(Cis) & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Trans-Interchromosomal',
    is.na(Cis) & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Trans-Intrachromosomal',
  )) 

co_eq_input <- co_eq_input %>% 
  mutate(
    log_exp = log(expression),
    log_exp = case_when(
      is.infinite(log_exp) ~ 0,
      TRUE ~ as.numeric(log_exp)
    )
  )
co_eq_input <-  distinct(co_eq_input)
p <- ggplot(data=co_eq_input, aes(x=interaction, y=log_exp,  colour=interaction, fill=interaction)) + 
  geom_violin() + geom_boxplot(width=0.16, color="grey30") +
  labs(x="Interaction", y="(Log) Median Expression (TPM)")+ 
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling2")) + 
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling2")) +
  theme_classic()+
  theme(legend.position = 'bottom', text = element_text(size = 14))

my_comparisons <- list(c('Cis', 'Trans-Intrachromosomal'))

p + facet_wrap(vars(tissue)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 10, size = 2.9) +stat_summary(fun=mean, geom="point", shape=17, size=1.2, colour='firebrick2')+
  geom_text(data = . %>% count(interaction, tissue), aes(y=-4.1, x=interaction,label=n), colour='black')  






#violin plot of cis eQTL expression, grouped by specific codes3d, shared, and specific GTEX
#isolating cis genes from codes3d
cis_co_input <- as_tibble(read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE)) %>% 
  dplyr::select(gencode_id, expression, tissue, interaction_type)
filter(cis_co_input, interaction_type == "Cis")
cis_co_input <- cis_co_input %>% dplyr::select(gencode_id, tissue, expression)
distinct(cis_co_input)
#isolating cis genes from GTEX with the right tissues
gtex_input <- as_tibble(read.csv(file='gtex_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz', sep='\t', header = TRUE, skip = 2))
gtex_results <- as_tibble(read.csv('results/gtex_results.txt', sep='\t', header=TRUE))

gtex_input <- gtex_input %>% dplyr::select(Name, Artery...Aorta, Artery...Coronary, Artery...Tibial, Heart...Atrial.Appendage, Heart...Left.Ventricle, Whole.Blood) %>% 
   pivot_longer(!Name, names_to="tissue", values_to="expression") %>% 
  mutate(tissue= case_when(
    tissue == "Artery...Aorta" ~ "Artery_Aorta",
    tissue == "Artery...Coronary" ~ "Artery_Coronary",
    tissue == "Artery...Tibial" ~ "Artery_Tibial",
    tissue == "Heart...Atrial.Appendage" ~ "Heart_Atrial_Appendage",
    tissue == "Heart...Left.Ventricle" ~ "Heart_Left_Ventricle",
    tissue == "Whole.Blood" ~ "Whole_Blood"
  ))
full_gtex <- gtex_input
full_gtex <- full_gtex %>% mutate_at(vars(expression),~log(.)) %>% rename(gencode_id=Name) %>% filter_at(vars("expression"), all_vars(!is.infinite(.)))
full_gtex$all_genes_gtex=1


gtex_input <- merge(gtex_results, gtex_input, by.x=c("gene_id", "tissue"), by.y=c("Name", "tissue"), how='LEFT') %>% dplyr::select(gene_id, tissue, expression)

cis_co_input <- cis_co_input %>% mutate_at(vars(expression),~log(.))
cis_co_input$codes3d=1
gtex_input <- gtex_input %>% mutate_at(vars(expression),~log(.)) %>% rename(gencode_id=gene_id)
gtex_input$gtex=1 


violin_input <- full_join(cis_co_input, gtex_input)
violin_input <-  full_join(violin_input, full_gtex)
violin_input <- violin_input %>% mutate(group= case_when(
  codes3d == 1 & is.na(gtex) ~ "CoDeS3D Specific",
  codes3d == 1 & gtex == 1 ~ "Shared",
  is.na(codes3d) & gtex == 1 ~ "GTEX Specific",
  all_genes_gtex == 1 ~ "GTEX (all genes)"
))



#merge the codes3d input with the gtex input, violin grouped by distinct codes3d, distinct gtex, and shared eqtls
violin_input <- violin_input %>% 
  mutate(
    expression = case_when(
      is.infinite(expression) ~ 0,
      TRUE ~ as.numeric(expression)
    )
  ) 

violin_input$group <- factor(violin_input$group, levels=c('CoDeS3D Specific', 'Shared', 'GTEX Specific', "GTEX (all genes)"))
violin_input <- distinct(violin_input)
p <- ggplot(data=violin_input, aes(x=group, y=expression,  colour=group, fill=group)) + 
  geom_violin() + geom_boxplot(width=0.16, color="white") +
  labs(x="Groups", y="(Log) Median Expression (TPM)")+ 
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1")) + 
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest1")) +
  theme_classic() +
  theme(legend.position = 'bottom', text = element_text(size = 14)) 

my_comparisons <- list(c('CoDeS3D Specific','Shared'),
                       c('Shared','GTEX Specific'), c('CoDeS3D Specific', 'GTEX Specific'),
                       c('CoDeS3D Specific', 'GTEX (all genes)'), c('Shared','GTEX (all genes)'), c('GTEX (all genes)','GTEX Specific'))

p + facet_wrap(vars(tissue)) + #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+  
  stat_compare_means(label.y = 15)+ 
  stat_summary(fun=mean, geom="point", shape=17, size=1.2, colour='firebrick2')+
  geom_text(data = . %>% count(group, tissue), aes(y=-7.7, x=group,label=n),colour='black')

#draft <- violin_input %>% filter(tissue == 'Artery_Tibial' & group == 'CoDeS3D Specific')
#draft <- left_join(draft, read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE), by = 'gencode_id')
#draft1 <- draft %>% dplyr::select(gene) %>% distinct()




#Extract the LOEUF scores for all the genes output from CoDeS3D, 
#using the "pLoF Metrics by Gene TSV" file under "Constraint" in gnomAD

loeuf_input <- read.csv('gnomad.v2.1.1.lof_metrics.by_gene.txt',sep='\t',header=TRUE)
gene_input <- read.csv('codes3d00/significant_eqtls.txt', sep='\t',header=TRUE) %>% dplyr::select(gencode_id, interaction_type)
gene_input <-  distinct(gene_input) %>% mutate(gene_id=gsub("\\..*","",gencode_id))
gene_input$present=1

gene_input <- gene_input %>% pivot_wider(values_from = present, names_from = interaction_type) %>%
  mutate(interaction = case_when(
    Cis == 1 & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans',
    Cis == 1 & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Cis & Trans-Interchromosomal', 
    Cis == 1 & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans-Intrachromosomal',
    is.na(Cis) & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Trans',
    Cis == 1 & is.na(`Trans-interchromosomal`) & is.na(`Trans-intrachromosomal`) ~ 'Cis',
    is.na(Cis) & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Trans-Interchromosomal',
    is.na(Cis) & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Trans-Intrachromosomal',
  )) 

co_loeuf_input <- merge(gene_input, loeuf_input, by.x="gene_id", by.y="gene_id", how='INNER')
co_loeuf_input <- distinct(co_loeuf_input)
p <- ggplot(data=co_loeuf_input, aes(x=interaction, y=oe_lof_upper,  colour=interaction, fill=interaction)) + 
  geom_violin() + geom_boxplot(width=0.1, color="white") +
  labs(title="CoDeS3D", x="Interaction", y="LOEUF Score")+
  scale_color_manual(values=wes_palette(n=6, name="BottleRocket1")) +
  scale_fill_manual(values=wes_palette(n=6, name="BottleRocket1")) +
  geom_text(data = . %>% count(interaction), aes(y=-0.2, x=interaction,label=n), colour='black')+
  theme(legend.position = 'bottom', text = element_text(size = 20))

my_comparisons <- list(c('Trans-Interchromosomal', 'Trans-Intrachromosomal'),
                       c('Cis & Trans-Intrachromosomal', 'Trans-Interchromosomal'),
                       c('Cis', 'Trans-Interchromosomal'))
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 3) +stat_summary(fun=mean, geom="point", shape=23, size=2, colour='white')


#LOEUF scores of http://130.216.236.202:8282/graphics/plot_zoom_png?width=1400&height=832all genes, codes3d, gtex, and shared
#cis genes from codes3d

cis_co_input_loeuf <- cis_co_input %>% mutate(gene_id=gsub("\\..*","",gencode_id)) %>% dplyr::select(gene_id, codes3d)
distinct(cis_co_input_loeuf)
gtex_input_loeuf <- gtex_input %>% mutate(gene_id=gsub("\\..*","",gencode_id)) %>% dplyr::select(gene_id, gtex)
distinct(gtex_input_loeuf)
all_gtex_loeuf <- full_gtex %>% mutate(gene_id=gsub("\\..*","",gencode_id)) %>% dplyr::select(gene_id, all_genes_gtex)
distinct(all_gtex_loeuf)

methods_loeuf_input <- full_join(cis_co_input_loeuf, gtex_input_loeuf) 
methods_loeuf_input <- full_join(methods_loeuf_input, all_gtex_loeuf)
methods_loeuf_input <- methods_loeuf_input %>% mutate(group= case_when(
  codes3d == 1 & is.na(gtex) ~ "CoDeS3D Specific",
  codes3d == 1 & gtex == 1 ~ "Shared",
  is.na(codes3d) & gtex == 1 ~ "GTEX Specific",
  all_genes_gtex == 1 ~ "GTEX (all genes)"
))
methods_loeuf_input <- distinct(methods_loeuf_input) %>% dplyr::select(gene_id, group)

methods_loeuf_input <- merge(methods_loeuf_input, loeuf_input, by.x='gene_id', by.y='gene_id', how='INNER')

methods_loeuf_input$group <- factor(methods_loeuf_input$group, levels=c('CoDeS3D Specific', 'Shared', 'GTEX Specific', 'GTEX (all genes)'))
p <- ggplot(data=methods_loeuf_input, aes(x=group, y=oe_lof_upper, colour=group, fill=group)) +
  geom_violin(trim=TRUE)+ geom_boxplot(width=0.1, colour='white') +
  labs(title="Cis eQTLs", x="Groups", y="LOEUF Score")+ 
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1")) + 
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest1")) +
  geom_text(data = . %>% count(group), aes(y=0, x=group,label=n),colour='black')+
  theme(legend.position = 'bottom', text = element_text(size = 20))

my_comparisons <- list(c('CoDeS3D Specific','Shared'),
                      c('CoDeS3D Specific', 'GTEX Specific'),
                       c('CoDeS3D Specific', 'GTEX (all genes)'), c('Shared','GTEX (all genes)'))

p + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+  
  stat_compare_means(label.y = 3) +stat_summary(fun=mean, geom="point", shape=23, size=2, colour='white')



#tissue specific LOEUF score violin plots
#codes3d eqtls
gene_input <- read.csv('codes3d00/significant_eqtls.txt', sep='\t',header=TRUE) %>% dplyr::select(gencode_id, interaction_type, tissue)
gene_input <-  distinct(gene_input, gencode_id, tissue, interaction_type) %>% mutate(gene_id=gsub("\\..*","",gencode_id))
gene_input$present=1

gene_input <- gene_input %>% pivot_wider(values_from = present, names_from = interaction_type) %>%
  mutate(Interaction = case_when(
    Cis == 1 & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans',
    Cis == 1 & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Cis & Trans-Inter', 
    Cis == 1 & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans-Intra',
    is.na(Cis) & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Trans',
    Cis == 1 & is.na(`Trans-interchromosomal`) & is.na(`Trans-intrachromosomal`) ~ 'Cis',
    is.na(Cis) & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Trans-Inter',
    is.na(Cis) & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Trans-Intra',
  )) 

co_loeuf_input <- merge(gene_input, loeuf_input, by.x="gene_id", by.y="gene_id", how='INNER')
co_loeuf_input <- distinct(co_loeuf_input)
p <- ggplot(data=co_loeuf_input, aes(x=Interaction, y=oe_lof_upper,  colour=Interaction, fill=Interaction)) + 
  geom_violin() + geom_boxplot(width=0.16, color="grey30") +
  labs(x="Interaction", y="LOEUF Score")+ 
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling2")) + 
  scale_fill_manual(values=wes_palette(n=4, name="Darjeeling2")) +
  theme_classic()+
  theme(legend.position = 'bottom', text = element_text(size = 14))

#my_comparisons <- list(c('Trans-Interchromosomal', 'Trans-Intrachromosomal'),
#                       c('Cis & Trans-Intrachromosomal', 'Trans-Interchromosomal'),c('Cis & Trans-Intrachromosomal', 'Trans-Intrachromosomal'),
#                       c('Cis', 'Trans-Intrachromosomal'), c('Cis', 'Cis & Trans-Intrachromosomal'),c('Cis', 'Trans-Interchromosomal'))

p + facet_wrap(vars(tissue)) + 
#  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 2.5) +stat_summary(fun=mean, geom="point", shape=17, size=1.2, colour='firebrick2')+
  geom_text(data = . %>% count(Interaction, tissue), aes(y=-0.1, x=Interaction,label=n), colour='black')



#tissue specific LOEUF score for cis eQTLs
gene_input <-  distinct(violin_input, gencode_id, tissue, group) %>% mutate(gene_id=gsub("\\..*","",gencode_id))
cis_loeuf_input <- merge(gene_input, loeuf_input, by.x="gene_id", by.y="gene_id", how='INNER')
cis_loeuf_input <- distinct(cis_loeuf_input)
p <- ggplot(data=cis_loeuf_input, aes(x=group, y=oe_lof_upper,  colour=group, fill=group)) + 
  geom_violin() + geom_boxplot(width=0.16, color="white") +
  labs(x="group", y="LOEUF Score")+ 
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1")) + 
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest1")) +
  theme_classic()+
  theme(legend.position = 'bottom', text = element_text(size = 14))

#my_comparisons <- list(c('CoDeS3D Specific','Shared'),
#                       c('Shared','GTEX Specific'), c('CoDeS3D Specific', 'GTEX Specific'),
#                       c('CoDeS3D Specific', 'GTEX (all genes)'), c('Shared','GTEX (all genes)'), c('GTEX (all genes)','GTEX Specific'))

p + facet_wrap(vars(tissue)) + 
#  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 2.5) +
  stat_summary(fun=mean, geom="point", shape=17, size=1.2, colour='firebrick2')+
  geom_text(data = . %>% count(group, tissue), aes(y=-0.1, x=group,label=n), colour='black')

