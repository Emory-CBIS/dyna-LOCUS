#!/bin/bash

# when you execute this script on the cluster run this command
# chmod u+x merge_decomp.sh
# ./merge_decomp.sh merge_decomp.R merge_decomp

##################### Change these constants ##############################
analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="/home/jran2/dFC/simulation/merge_decomp/scratch"  # location of your persistent scratch dir
resultdir="/home/jran2/dFC/simulation/merge_decomp/scratch/out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
max_jobs=42   # max number of jobs to run at a time
total_jobs=42   # total number of jobs
############## typically you don't have to change anything below here #######

username=$(id -nu)

nloops=$((${total_jobs}/${max_jobs}-1))

# run R script to get list size
# this apparently doesn't work
# listsize=$(${script} listsize)

# if scratch directory doesn't exist, make it
[ ! -d ${myscratch} ] && mkdir ${myscratch}
[ ! -d ${myscratch}/out ] && mkdir ${myscratch}/out
[ ! -d ${myscratch}/err ] && mkdir ${myscratch}/err

# submit first batch of jobs
for i in $(seq 1 ${max_jobs}); do
	echo "#!/bin/bash" >> script$i.sh
#	echo "#SBATCH --nodes=1 # ask for 1 node" >> script$i.sh
#	echo "#SBATCH --ntasks-per-node=4 # 4 tasks each node" >> script$i.sh
#	echo "#SBATCH --mem-per-cpu=10G" >> script$i.sh
	echo "#SBATCH --partition=preemptable" >> script$i.sh
	echo "#SBATCH --job-name=${analysis}$i" >> script$i.sh
	echo "#SBATCH --error=${myscratch}/err/${analysis}$i.err" >> script$i.sh
	echo "#SBATCH --output=${myscratch}/out/${analysis}$i.out" >> script$i.sh

    echo "source ~/miniconda3/etc/profile.d/conda.sh" >> script$i.sh
	echo "conda activate rproj" >> script$i.sh

	echo "Rscript ${script} $i" >> script$i.sh
    
	if [[ $i -ne ${max_jobs} ]]; then
	   sbatch script$i.sh
	else
	   holdid=$(sbatch script$i.sh | sed 's/Submitted batch job //')
    fi
done

# sub next batches of jobs holding 
for j in $(seq 1 $nloops); do
	for i in $(seq 1 ${max_jobs}); do
		# hold_name=${analysis}$((($j-1)*${max_jobs}+$i))
		jid=$(($j*${max_jobs}+$i))

		echo "#!/bin/bash" >> script$jid.sh
#		echo "#SBATCH --nodes=1 # ask for 1 node" >> script$jid.sh
#		echo "#SBATCH --ntasks-per-node=4 # 4 tasks each node" >> script$jid.sh
#		echo "#SBATCH --mem-per-cpu=10G" >> script$jid.sh
		echo "#SBATCH --partition=preemptable" >> script$jid.sh
		echo "#SBATCH --job-name=${analysis}$jid" >> script$jid.sh
		echo "#SBATCH --error=${myscratch}/err/${analysis}$jid.err" >> script$jid.sh
		echo "#SBATCH --output=${myscratch}/out/${analysis}$jid.out" >> script$jid.sh

		echo "source ~/miniconda3/etc/profile.d/conda.sh" >> script$jid.sh
		echo "conda activate rproj" >> script$jid.sh

		echo "Rscript ${script} $jid" >> script$jid.sh

		if [[ $i -ne ${max_jobs} ]]; then
	   		sbatch --dependency=afterany:$holdid script$jid.sh
		else
	   		holdid=$(sbatch --dependency=afterany:$holdid script$jid.sh | sed 's/Submitted batch job //')
    	fi
	done
done