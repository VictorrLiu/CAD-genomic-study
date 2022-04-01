library(dplyr)
library(UpSetR)

#using upsetR to show the genes overlap between differnt tissues
upset_input <- read.csv('codes3d00/significant_eqtls.txt', sep='\t', header = TRUE)
genes_upset_input <- upset_input %>% 
  distinct(gene, tissue) %>% 
  mutate(present=1) %>% 
  pivot_wider(names_from = tissue, values_from = present) %>% 
  dplyr::select(!gene) 

genes_upset_input[is.na(genes_upset_input)]=0
genes_upset_input %>% 
   as.data.frame() %>% 
   upset(order.by = 'freq', 
 mainbar.y.label = 'Number of Genes',
 sets.x.label = 'Genes Per Tissue', nintersects = NA,nsets =6,text.scale = c(3, 3, 2, 3, 3, 2.2), line.size = 1.5, point.size = 5) 


#same as above, but only for the arterial tissues
genes_artery_upset_input <- genes_upset_input %>% filter(Artery_Aorta == 1 | Artery_Coronary == 1 | Artery_Tibial == 1) %>% 
  dplyr::select('Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial')

upset(as.data.frame(genes_artery_upset_input),
      mainbar.y.label = 'Number of Genes',
      sets = c('Artery_Tibial', 'Artery_Coronary', 'Artery_Aorta'),keep.order = TRUE,
      sets.x.label = 'Genes Per Tissue', nintersects = NA, text.scale = c(1.5, 1.3, 1.2, 1, 1.5, 1.5), line.size = 0.8, point.size = 1.8)


#upsetR plot for eQTL SNPs overlap
#all six tissues
eqtls_upset_input <- upset_input %>% 
   distinct(snp, tissue) %>% 
   mutate(present=1) %>% 
   pivot_wider(names_from = tissue, values_from = present) %>% 
   dplyr::select(!snp) 

eqtls_upset_input[is.na(eqtls_upset_input)]=0
eqtls_upset_input %>% 
   as.data.frame() %>% 
   upset(order.by = 'freq', 
         mainbar.y.label = 'Number of eQTL SNPs',
         sets.x.label = 'eQTL SNPs Per Tissue', nintersects = NA,nsets =6, text.scale = c(3, 3, 2, 3, 3, 2.2), line.size = 1.5, point.size = 5) 


#same as above, but only for the arterial tissues
eqtls_artery_upset_input <- eqtls_upset_input %>% filter(Artery_Aorta == 1 | Artery_Coronary == 1 | Artery_Tibial == 1) %>% 
   dplyr::select('Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial')

upset(as.data.frame(eqtls_artery_upset_input), order.by = 'freq', 
      mainbar.y.label = 'Number of eQTL', 
      sets = c('Artery_Tibial', 'Artery_Coronary', 'Artery_Aorta'),keep.order = TRUE,
      sets.x.label = 'SNPs Per Tissue', nintersects = NA, text.scale = c(1.5, 1.3, 1.2, 1, 1.5, 1.5), line.size = 0.8, point.size = 1.8)


#upsetR plot for eQTL Gene Pairs overlap
#all six tissues
pairs_upset_input <- upset_input %>% 
   distinct(snp, gene, tissue) %>% 
   mutate(present=1) %>% 
   pivot_wider(names_from = tissue, values_from = present) #%>% 
  # dplyr::select(-c(snp, gene)) 

pairs_upset_input[is.na(pairs_upset_input)]=0
pairs_upset_input %>% 
   as.data.frame() %>% 
   upset(order.by = 'freq', 
         mainbar.y.label = 'Number of eQTL Gene Pairs',
         sets.x.label = 'eQTL Gene Pairs Per Tissue', nintersects = NA,nsets =6, text.scale = c(3, 3, 2, 3, 3, 2.2), line.size = 1.5, point.size = 5) 


#same as above, but only for the arterial tissues
pairs_artery_upset_input <- pairs_upset_input %>% filter(Artery_Aorta == 1 | Artery_Coronary == 1 | Artery_Tibial == 1) %>% 
   dplyr::select('Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial')

upset(as.data.frame(pairs_artery_upset_input), order.by = 'freq', 
      mainbar.y.label = 'Number of eQTL',
      sets.x.label = 'Pairs Per Tissue', 
      sets = c('Artery_Tibial', 'Artery_Coronary', 'Artery_Aorta'),keep.order = TRUE,
      nintersects = NA, text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.5), line.size = 0.8, point.size = 1.8)


#coronary and tibial pairs that are overlapped.

overlap_examples <- upset_input %>% 
   distinct(snp, gene, tissue) %>% 
   mutate(present=1) %>% 
   pivot_wider(names_from = tissue, values_from = present) %>% 
   filter(Artery_Coronary == 1 & Artery_Tibial == 1 & is.na(Artery_Aorta))

print(overlap_examples)

