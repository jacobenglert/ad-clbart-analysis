#!/bin/bash
#SBATCH --job-name=ad-hw
#SBATCH --nodes=1                       # Node count (number of computers)
#SBATCH --ntasks-per-node=1             # total number of tasks across all nodes
#SBATCH --cpus-per-task=1               # cpu-cores per task (>1 if multi-threaded tasks, internal parallelization within R for example)
#SBATCH --mem-per-cpu=6G                # memory per cpu-core (4G is default)
#SBATCH --time=7-00:00:00               # total run time limit (HH:MM:SS)
#SBATCH --array=1-200             	 # Job array with index values. Show up as SLURM_ARRAY_TASK_ID
#SBATCH --partition=chang        # Designation of which cluster nodes to launch jobs on
#SBATCH --mail-type=FAIL                # Email received when job fails only (can also set as start, etc.)
#SBATCH --mail-user=jrengle@emory.edu   # Email to receive
#SBATCH --output=Output/Heatwave/slurm-%A.%a.out # Output file
#SBATCH --error=Error/Heatwave/slurm-%A.%a.err   # Error file

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID."
echo "Executing on the machine:" $(hostname)

module load R/4.2.2

Rscript R/ad-hw-fit.R $SLURM_ARRAY_TASK_ID
