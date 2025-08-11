#!/bin/bash

#SBATCH --job-name=ALAN_T
#SBATCH --output=/work/%u/%x-%A-%a.out
#SBATCH --time=0-10:00:00
#SBATCH --mem-per-cpu=2G

module load foss/2022b R/4.2.2

output_dir=/work/$USER/ALAN_T
mkdir -p $output_dir

cd /home/$USER/ALAN_Temp_foodwebs/code

Rscript \
  --vanilla \
  HPC_ALAN_T_code.R \
  $params \
  $output_dir