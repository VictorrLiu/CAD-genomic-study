#!/usr/bin/bash
in_dir="/mnt/projects/users/vliu378/project2/results/final_results/coronary_genes_input.txt"
out_dir="/mnt/projects/users/vliu378/project2/results/final_results/coronary_output"
mkdir "/mnt/projects/users/vliu378/project2/results/final_results/coronary_output"
python3 /mnt/projects/ecomorbidity/query_string.py \
	-i $in_dir \
	-o $out_dir/coronary_string_ppin.txt \
	-l 4
python3 /mnt/projects/ecomorbidity/get_ppi_eqtls.py \
	-p $out_dir/coronary_string_ppin.txt \
	-o $out_dir \
	-e /mnt/projects/tissue_maps/artery_coronary/Artery_Coronary/
python3 /mnt/projects/ecomorbidity/find_snp_disease.py \
	-e $out_dir \
	-o $out_dir	
