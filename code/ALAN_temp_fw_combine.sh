#!/bin/bash

#SBATCH --job-name=HPC_ALAN_temp_fw_combine
#SBATCH --output=/work/%u/%u-%A-%a.out
#SBATCH --time=0-00:05:00
#SBATCH --mem-per-cpu=1G

module load foss/2022b R/4.2.2

cd /home/$USER/ALAN_Temp_foodwebs/code

Rscript \
  --vanilla \
  HPC_ALAN_temp_fw_combine.R
