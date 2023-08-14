#!/bin/bash

# when you execute this script on the cluster run this command
# chmod u+x decomposition_DL.sh
# ./decomposition_DL.sh decomposition_DL.py decomp_DL

##################### Change these constants ##############################
analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="/home/jran2/dFC/simulation/decomp_DL/scratch"  # location of your persistent scratch dir
resultdir="/home/jran2/dFC/simulation/decomp_DL/scratch/out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
total_jobs=3   # total number of jobs
############## typically you don't have to change anything below here #######

username=$(id -nu)

# if scratch directory doesn't exist, make it
[ ! -d ${myscratch} ] && mkdir ${myscratch}
[ ! -d ${myscratch}/out ] && mkdir ${myscratch}/out
[ ! -d ${myscratch}/err ] && mkdir ${myscratch}/err

# submit first batch of jobs
for i in $(seq 1 ${total_jobs}); do
	echo "#!/bin/bash" >> script$i.sh
#	echo "#SBATCH --nodes=1 # ask for 1 node" >> script$i.sh
#	echo "#SBATCH --ntasks-per-node=4 # 4 tasks each node" >> script$i.sh
#	echo "#SBATCH --mem-per-cpu=10G" >> script$i.sh
	echo "#SBATCH --partition=preemptable" >> script$i.sh
	echo "#SBATCH --job-name=${analysis}$i" >> script$i.sh
	echo "#SBATCH --error=${myscratch}/err/${analysis}$i.err" >> script$i.sh
	echo "#SBATCH --output=${myscratch}/out/${analysis}$i.out" >> script$i.sh

    echo "iter=$i" >> script$i.sh
	echo "export iter" >> script$i.sh

	echo "source ~/miniconda3/etc/profile.d/conda.sh" >> script$i.sh
	echo "conda activate rproj" >> script$i.sh

	echo "python ${script}" >> script$i.sh
    
	if [[ $i -ne ${max_jobs} ]]; then
	   sbatch script$i.sh
	else
	   holdid=$(sbatch script$i.sh | sed 's/Submitted batch job //')
    fi
done