import pandas as pd
import numpy as np 

total_genes = pd.read_csv('/mnt/projects/users/vliu378/project1/codes3d00/significant_eqtls.txt', sep='\t')
df = pd.DataFrame(total_genes, columns = ['snp', 'gene', 'tissue'])

df_list = df.loc[(df['tissue'] == 'Artery_Coronary') | (df['tissue'] == 'Artery_Tibial')]

output_list = df_list[df_list.duplicated(['snp', 'gene'])]
keep_last = df_list[df_list.duplicated(['snp', 'gene'], keep='last')]
output_list.append(keep_last)
print(keep_last)

#output_list.to_csv('/mnt/projects/users/vliu378/project1/results/c.txt', sep= '\t', index = False) 
