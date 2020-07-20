var_array=( "0.01" "0.02" "0.03" "0.04" "0.05" )
pi_array=( "0.0" "0.25" "0.5" "0.75" "1.0" )
for i in "${var_array[@]}"
do
	for j in "${pi_array[@]}"
	do
		mkdir "var_ratio_${i}_pi_ratio_${j}"
	done
done
