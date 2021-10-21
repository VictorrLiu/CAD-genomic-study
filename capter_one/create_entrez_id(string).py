#Louvain analysis on the STRING analysis of the genes output from CoDeS3D (grouped by tissues)

#convert the gene names to Entrez_id for Louvain (myGene package)
import pandas as pd
import numpy as np
import mygene
mg = mygene.MyGeneInfo()

total_genes = pd.read_csv('/mnt/projects/users/vliu378/project1/codes3d00/significant_eqtls.txt', sep='\t')
total_genes = total_genes[['gene','tissue']].drop_duplicates()

out = mg.querymany(total_genes['gene'], scopes='symbol',fields='entrezgene',species='human')


entrez_id = pd.DataFrame.from_dict(out)
entrez_id.dropna(subset = ['entrezgene'], inplace=True)
final_total_genes = pd.merge(total_genes, entrez_id, how='inner', left_on='gene', right_on='query')
final_total_genes = final_total_genes[['gene', 'tissue', '_id']].drop_duplicates()

Artery_Aorta_id = final_total_genes[final_total_genes['tissue'] == 'Artery_Aorta']
Artery_Coronary_id = final_total_genes[final_total_genes['tissue'] == 'Artery_Coronary']
Artery_Tibial_id = final_total_genes[final_total_genes['tissue'] == 'Artery_Tibial']
Heart_Atrial_Appendage_id = final_total_genes[final_total_genes['tissue'] == 'Heart_Atrial_Appendage']
Heart_Left_Ventricle_id = final_total_genes[final_total_genes['tissue'] == 'Heart_Left_Ventricle']
Whole_Blood_id = final_total_genes[final_total_genes['tissue'] == 'Whole_Blood']

Artery_Aorta_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Artery_Aorta_id.txt', index=False, header=False)
Artery_Coronary_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Artery_Coronary_id.txt', index=False, header=False)
Artery_Tibial_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Artery_Tibial_id.txt', index=False, header=False)
Heart_Atrial_Appendage_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Heart_Atrial_Appendage_id.txt', index=False, header=False)
Heart_Left_Ventricle_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Heart_Left_Ventricle_id.txt', index=False, header=False)
Whole_Blood_id[['_id']].to_csv('/mnt/projects/users/vliu378/project1/script/entrez_id/Whole_Blood_id.txt', index=False, header=False)