#!/bin/bash
#SBATCH -A HENDERSON-SL2-CPU
#SBATCH -p icelake
#SBATCH --job-name=pore_c_BglII_TAIR12
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00

# Load required modules
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl

# Define input and output base directories


bash /rds/project/rds-CXkCVqbkvMY/KJ_WORKING/2025April/run_porec_analysis.sh $1 $2


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to pwd.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: date"
echo "Running on master node: hostname"
echo "Current directory: pwd"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=generate_pbs_nodefile
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo cat machine.file.$JOBID | sed -e 's/\..*$//g'
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
