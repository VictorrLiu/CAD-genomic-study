#!/usr/bin/bash
sim_dir=/mnt/projects/users/vliu378/project2/results/coronary/coronary_bootstraps_simulations
out_dir=/mnt/projects/users/vliu378/project2/results/coronary/coronary_bootstraps_enrichment
for sim in $(ls $sim_dir)
do
	i=${sim::-4}
	echo "simulation $i"
	mkdir -p $out_dir/$i
	python3 /mnt/projects/ecomorbidity/query_string.py \
		-i $sim_dir/$sim \
	       	-o $out_dir/$i/string_ppin.txt \
		-l 4
	python3 /mnt/projects/ecomorbidity/get_ppi_eqtls.py \
		-p $out_dir/$i/string_ppin.txt \
		-o $out_dir/$i/ \
		-e /mnt/projects/tissue_maps/artery_coronary/Artery_Coronary/
	python3 /mnt/projects/ecomorbidity/find_snp_disease.py \
		-e $out_dir/$i/ \
		-o $out_dir/$i/	
done
