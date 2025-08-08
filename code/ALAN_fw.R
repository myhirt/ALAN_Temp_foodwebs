library(ATNr)
library(tidyverse)
library(parallel)
library(rstudioapi)
library(purrr)
library(furrr)

setwd(dirname(getActiveDocumentContext()$path))

source("ALAN_fw_functions.R")

# set basic parameters for simulations
n_species <- 60
n_basal <- 20
n_cons = n_species - n_basal
n_nut <- 2
diel_groups = c('D', 'N', 'C') # possible diel groups for species
probs = c(1,1,1) # distribution of the different diel groups across the community

n_reps <- 100

# set community parameters (body masses and diel groups
community_list <- replicate(n_reps, {
  masses <- 10 ^ c(
    sort(runif(n_basal, 0, 3)),
    sort(runif(n_cons, 2, 5))
  )
  biomasses <- runif(n_species + n_nut, 2, 3)
  L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)
  fw <- L
  fw[fw > 0] <- 1
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  diel_group_fw <- sample(diel_groups, n_cons, replace = TRUE, prob = probs)

  list(model = model,
       diel_group_fw = diel_group_fw)
}, simplify = FALSE)

# Define parameter combinations
param_grid <- expand.grid(
  light.scenario = c("N", "C", "C.N", "C.D"),
  index = 1:n_reps,
  temp = c(10, 15),
  x = seq(0, 0.5, by = 0.1),
  stringsAsFactors = FALSE
)

# Attach corresponding community to each row based on `index`
param_grid$community <- community_list[param_grid$index]

# parallelize
plan(multisession, workers = 30) 

# Run the simulation for each parameter combination
results <- future_pmap_dfr(param_grid, function(light.scenario, index, temp, x, community) {
  df <- sim.ALAN(
    x = x,
    model = community$model,
    light.scenario = light.scenario,
    diel_group = community$diel_group_fw,
    temp = temp
  )
  df$index <- index
  df$light.scenario <- light.scenario
  return(df)
  },  
.options = furrr_options(seed = TRUE)
)

ggplot(results, aes(x=x, y = tot.pers, color = light.scenario)) +
  geom_jitter()+
  geom_smooth(method = "loess") +
  labs(x = "light pollution") +
  theme_minimal()

ggplot(results, aes(x=x, y = tot.pers, color = light.scenario)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~temp) +
  labs(x = "light pollution") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file = "output/scenario_all.pdf", height=6, width =8)

ggplot(filter(results, light.scenario == "N"), aes(x=x, y = tot.pers)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~temp) +
  labs(x = "light pollution") +
  ggtitle("scenario nocturnal") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file = "../output/scenario_N.pdf", height=6, width =8)

ggplot(filter(results, light.scenario == "C"), aes(x=x, y = tot.pers)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~temp) +
  labs(x = "light pollution") +
  ggtitle("scenario crepuscular") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file = "../output/scenario_C.pdf", height=6, width =8)


ggplot(filter(results, light.scenario == "C.N"), aes(x=x, y = tot.pers)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~temp) +
  labs(x = "light pollution") +
  ggtitle("scenario crepuscular & nocturnal") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file = "../output/scenario_C.N.pdf", height=6, width =8)

ggplot(filter(results, light.scenario == "N.D"), aes(x=x, y = tot.pers)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~temp) +
  labs(x = "light pollution") +
  ggtitle("scenario diurnal & nocturnal") +
  theme_minimal() + theme(legend.position = "none")

ggsave(file = "../output/scenario_N.D.pdf", height=6, width =8)



######
results.old <- read.csv("../data/results_new.csv")

ggplot(filter(results.old, temperature == 15), aes(x=light, y = tot_ext, color = light.effect)) +
  #geom_point()+
  geom_smooth(method = "lm") +
  labs(x = "light pollution") +
  theme_minimal() +
  facet_wrap(~S)

unique(results.old$temperature)
