#Inner-merge (input)"snps.txt" and Gtex files of different tissues, so only the lines with matching variant_id are extracted

import numpy as np
import pandas as pd
import os
import sys

snps_fp = '/mnt/projects/users/vliu378/project1/codes3d00/snps.txt'
tissues_fp = '/mnt/projects/users/vliu378/project1/gtex_data'
output_fp = '/mnt/projects/users/vliu378/project1/results/gtex_results.txt'

tissues = ['Artery_Aorta', 'Artery_Coronary','Artery_Tibial', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Whole_Blood']

snps = pd.read_csv(snps_fp, sep='\t')
snps= snps[['snp','variant_id']].drop_duplicates()

results = []
for tissue in tissues:
  print(tissue)
  df = pd.read_csv(os.path.join(tissues_fp, tissue + '.v8.signif_variant_gene_pairs.txt.gz'), sep = '\t', compression = 'gzip')
  df['tissue'] = tissue
  results.append(snps.merge(df, how = "inner", on = "variant_id"))
    
results = pd.concat(results)

results.to_csv(output_fp, sep= '\t', index = False) 

#Same thing but to Trans-Gtex data
trans_output_fp = '/mnt/projects/users/vliu378/project1/results/trans_gtex_results.txt'

trans_df = pd.read_csv('/mnt/projects/users/vliu378/project1/gtex_data/GTEx_Analysis_v8_trans_eGenes_fdr05.txt', sep='\t')
trans_results = snps.merge(df, how = "inner", on = "variant_id")

trans_results.to_csv(trans_output_fp, sep= '\t', index = False) 

#Number of eqtls and genes from the output Gtex data from above
gtex_fp = '/mnt/projects/users/vliu378/project1/results/gtex_results.txt'

gtex = pd.read_csv(gtex_fp, sep='\t')
print(gtex[['snp']].drop_duplicates())
print(gtex[['gene_id']].drop_duplicates())

snps_fp = '/mnt/projects/users/vliu378/project1/codes3d00/snps.txt'
snps = pd.read_csv(snps_fp, sep='\t')
snps= snps[['snp','variant_id']].drop_duplicates()
print(snps)
