#!/bin/bash
wc -l /humgen/diabetes2/users/lthakur/lap_test/ldsc/raw/cp_gene_sets/*.txt | sed 's;/humgen/diabetes2/users/lthakur/lap_test/ldsc/raw/cp_gene_sets/;;' | sed 's;.txt;;' | awk '$2 != "total" && $1 >= 10 {print $2,"class","gene_set"; print $2,"parent","euro"; print $2,"gene_set_gene_set_list_file","@raw_dir/cp_gene_sets/"$2".txt"}' > /home/unix/flannick/lap/projects/ldsc/meta/cp_gene_sets.meta

#wc -l /humgen/diabetes2/users/lthakur/lap_test/ldsc/raw/cp_gene_sets/*.txt | head
