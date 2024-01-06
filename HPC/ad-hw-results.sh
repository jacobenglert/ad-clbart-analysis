#!/bin/bash

array_job0=$(sbatch --parsable HPC/ad-hw-results.slurm)

sbatch --depend=afterany:$array_job0 HPC/ad-hw-results-combine.slurm
