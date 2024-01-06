#!/bin/bash
#SBATCH --job-name=ad-hw-combine      # create a name for the job
#SBATCH --nodes=1                         # Node count (number of computers)
#SBATCH --ntasks-per-node=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                 # cpu-cores per task (>1 if multi-threaded tasks, internal parallelization within R for example)
#SBATCH --mem-per-cpu=20G                  # memory per cpu-core (4G is default)
#SBATCH --partition=chang               # Designation of which cluster nodes to launch jobs on
#SBATCH --mail-type=ALL                   # Email received when job fails only (can also set as start, etc.)
#SBATCH --mail-user=jrengle@emory.edu     # Email to receive
#SBATCH --output=Output/Heatwave/combine.out       # Output file
#SBATCH --error=Error/Heatwave/combine.err         # Error file

module load R/4.2.2

Rscript R/ad-hw-results-combine-hpc.R
