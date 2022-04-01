library(tidyverse)
library(ggplot2)
library(wesanderson)
library(dbplyr)
library(dplyr)
library(ggpubr)

#absolute fold change violin plots (interaction types)

afc_input <- as_tibble(read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE))
afc_input <- afc_input %>% dplyr::select(snp, gene, log2_aFC, interaction_type, tissue)
afc_input$present=1

afc_input <- afc_input %>% pivot_wider(values_from = present, names_from = interaction_type) %>%
  mutate(interaction = case_when(
    Cis == 1 & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans',
    Cis == 1 & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Cis & Trans-Interchromosomal', 
    Cis == 1 & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Cis & Trans-Intrachromosomal',
    is.na(Cis) & `Trans-interchromosomal` == 1 & `Trans-intrachromosomal` == 1 ~ 'Trans',
    Cis == 1 & is.na(`Trans-interchromosomal`) & is.na(`Trans-intrachromosomal`) ~ 'Cis',
    is.na(Cis) & `Trans-interchromosomal` == 1 & is.na(`Trans-intrachromosomal`) ~ 'Trans-Interchromosomal',
    is.na(Cis) & is.na(`Trans-interchromosomal`) & `Trans-intrachromosomal` == 1 ~ 'Trans-Intrachromosomal',
  )) 
my_colours <- wes_palette(n=4, name="Darjeeling2")[c(1,3,4)]

p <- ggplot(data=afc_input, aes(x=interaction, y=log2_aFC,  colour=interaction, fill=interaction)) + 
  geom_violin() + geom_boxplot(width=0.16, color="grey30") +
  labs(x="Interaction", y="Allelic Fold Change (Log2)")+ 
  scale_color_manual(values=my_colours) + 
  scale_fill_manual(values=my_colours) +
  theme_classic()+
  theme(legend.position = 'bottom', text = element_text(size = 14))

my_comparisons <- list(c('Cis', 'Trans-Interchromosomal'),c('Trans-Interchromosomal', 'Trans-Intrachromosomal'),c('Cis', 'Trans-Intrachromosomal'))

p + facet_wrap(vars(tissue)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 2.5, size = 2.9) +
  stat_summary(fun=mean, geom="point", shape=17, size=1.2, colour='firebrick2')+
  geom_text(data = . %>% count(interaction, tissue), aes(y=-2.1, x=interaction,label=n), colour='black')  



