#!/bin/bash

#SBATCH --job-name=ALAN_temp_fw_original
#SBATCH --output=/work/%u/%x-%A-%a.out
#SBATCH --time=0-00:10:00
#SBATCH --mem-per-cpu=1G

module load foss/2022b R/4.2.2

output_dir=/work/$USER/ALAN_temp_fw_original
mkdir -p $output_dir

params=/home/$USER/ALAN_Temp_foodwebs/code/HPC_ALAN_temp_fw_params.csv

cd /home/$USER/ALAN_Temp_foodwebs/code

Rscript \
  --vanilla \
  HPC_ALAN_code_benoit_original.R \
  $params \
  $output_dir
