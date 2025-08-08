#------------ Load all packages ------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ATNr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(docopt))

doc <- "usage: HPC_ALAN_code_benoit_original.R <params> <output_dir>"
opts <- docopt(doc)

## read parameter file
params <- fread(opts$params)
output_dir <- opts$output_dir

## try to get task id
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

rep <- params[task]$rep
light.effect <- params[task]$light.effect
t1 <- params[task]$t1
t2 <- params[task]$t2
S <- params[task]$S
n_species <- params[task]$n_species
n_basal <- params[task]$n_basal
n_nut <- params[task]$n_nut
n_cons <- params[task]$n_species - params[task]$n_basal
masses <- as.numeric(strsplit(params$masses[task], ";")[[1]])
biomasses <- as.numeric(strsplit(params$biomasses[task], ";")[[1]])
period <- as.numeric(strsplit(params$period[task], ";")[[1]])

source("HPC_ALAN_code_benoit_original_functions.R")

set.seed(12)

res = run.light.gradient(light.effect=light.effect, 
                         t1=t1, 
                         t2=t2, 
                         S=S, 
                         n_basal=n_basal, 
                         n_species=n_species, 
                         n_nuts=n_nuts, 
                         masses = masses, 
                         biomasses=biomasses, 
                         period =period,
                         rep = rep)


write.csv(
  res,
  file.path(
    opts$output_dir,
    paste0("light_effect_", light.effect, "_", rep, ".csv")
  )
)



