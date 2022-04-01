#generating 10k simulations of 166 genes for bootstrap used 

import pandas as pd

ref_input = "/mnt/projects/codes3d/lib/reference_files/genes/gene_reference.bed"
coronary_bootstrap_intermediate_output = "/mnt/projects/users/vliu378/project2/results/coronary/coronary_bootstraps"
aorta_bootstrap_intermediate_output = "/mnt/projects/users/vliu378/project2/results/aorta/aorta_bootstraps"

ref = pd.read_csv(ref_input, sep='\t')
ref.columns = ['chromosome', 'ID1', 'ID2', 'gene', 'genecode']
ref_gene_list = ref['gene']
ref_gene_list = ref_gene_list.drop_duplicates() 


for i in range(1, 10001):
  seed = i
  a = ref_gene_list.sample(n = 166, random_state = seed)
  a.to_csv((coronary_bootstrap_intermediate_output + "/" + str(seed) + ".txt"), sep='\t', index=False)

  b = ref_gene_list.sample(n = 346, random_state = seed)
  b.to_csv((aorta_bootstrap_intermediate_output + "/" + str(seed) + ".txt"), sep='\t', index=False)