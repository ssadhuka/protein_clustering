sleeper=0 
NUM_TO_SUBMIT=999
CAPACITY=1000
tac proteins.txt | while read line || [[ -n $line ]] ;
do
	#for i in {1..${NUM_TO_SUBMIT}}
	#do
	echo "source /broad/software/scripts/useuse" > job_scripts/submit_${line}.sh
	chmod u+w job_scripts/submit_${line}.sh
	echo "reuse UGER" >> job_scripts/submit_${line}.sh
	echo "reuse Anaconda" >> job_scripts/submit_${line}.sh	
	echo "python /humgen/diabetes2/users/ssadhuka/ssadhuka/protein_clusters/Scale_Up/uniprot_tagger_copy.py --uniprot ${line}" >> job_scripts/submit_${line}.sh
	#python /humgen/diabetes2/users/ssadhuka/ssadhuka/protein_clusters/Scale_Up/uniprot_tagger_copy.py --uniprot "${line}"
	#echo ${file}
	qsub job_scripts/submit_${line}.sh -N ${line} -o outputs_uniprot/${line}.out -j y
	# python uniprot_tagger_copy.py --uniprot "${line}"
	#done;
	NUM_RUNNING=$(qstat | wc -l)
	NUM_TO_SUBMIT="$((CAPACITY-NUM_RUNNING))"
	
	while (($NUM_TO_SUBMIT <= 100))
	do
		NUM_RUNNING=$(qstat | wc -l)
		NUM_TO_SUBMIT="$((CAPACITY-NUM_RUNNING))"
	done;
	
	#echo ${NUM_TO_SUBMIT}
	#if (( ${sleeper} % 1000 == 999 ))
	#then
	#	echo "sleeping"
	#	sleep 6m
	#fi
	#((sleeper++))
done;
