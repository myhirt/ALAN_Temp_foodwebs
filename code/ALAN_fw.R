library(ATNr)
library(tidyverse)
library(parallel)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
set.seed(12)

source("ALAN_fw_functions.R")

# set parameters for simulations
S.all = c(1, 2)
n_species <- 60
n_basal <- 20
n_cons = n_species - n_basal
n_nut <- 2
scenarios = c('N', 'C', 'C.N', 'N.D') # light pollution scenarios, which groups are affected
diel_groups = c('D', 'N', 'C') # possible diel groups
probs = c(1,1,1) 

masses <- 10 ^ c(sort(runif(n_basal, 0, 3)),
                 sort(runif(n_species - n_basal, 2, 5)))
biomasses <- runif(n_species + n_nut, 2, 3)
diel_group_fw = sample(diel_groups, n_cons, replace = TRUE, prob = probs)

light.scenario = "CN"

test.ALAN <- NULL  # initialize before loop

for (i in 1:10) {  # run 100 times
  for (temp in seq(0, 40, by = 5)) {
    for (x in seq(0, 0.5, by = 0.1)) {
      test.ALAN.new <- sim.ALAN(x, masses, n_basal, light.scenario, diel_group_fw, temp)  
      test.ALAN <- rbind(test.ALAN, test.ALAN.new)
    }
  }
}

ggplot(test.ALAN, aes(x=x, y = tot.ext, color = as.factor(temp))) +
  geom_point()+
  geom_smooth(method = "lm")


ggplot(test.ALAN, aes(x=temp, y = tot.ext, color = as.factor(x))) +
  geom_point()+
  geom_smooth(method = "loess")








run.light.gradient = function(param){
  # run a light intensity gradient at 2 temperatures, 
  
  # create the L matrix
  L <- create_Lmatrix(param$masses, param$n_basal, Ropt = 100, gamma = 2, th = 0.01)
  # create the 0/1 version of the food web
  fw <- L
  fw[fw > 0] <- 1
  model <- create_model_Unscaled_nuts(param$n_species, param$n_basal, param$n_nut, param$masses, fw)
  
  light = seq(0, 0.5, by = 0.05)
  model <- initialise_default_Unscaled_nuts(model, L, temperature = param$t1)
  model$S = rep(param$S, param$n_nut)
  model$q = rep(1.2, n_species - n_basal)
  
  # test if the separation btween nocturnal and diurnal lead to some new basal species
  # even without light scenario
  # if yes, cancel run
  overlaps = overlap.modification(x, diel_group_fw, 0)
  b2 = model$b
  b2[(n_basal +1):n_species, ] = b2[(n_basal +1):n_species, ] * overlaps
  if (any(colSums(b2) == 0.0)){return(NULL)}
  
  exts.t1 = sapply(light, run.light, model, 
                     light.scenario = param$light.scenario, diel_group_fw = param$diel_group_fw)
  res1 = cbind.data.frame(t(exts.t1), light, param$t1, param$light.scenario, param$S, param$rep)
  names(res1) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x", 
                  "light", "temperature", "light.scenario", "S", "replicate")
  
  model <- initialise_default_Unscaled_nuts(model, L, temperature = param$t2)
  model$S = rep(param$S, param$n_nut)
  model$q = rep(1.2, n_species - n_basal)
  exts.t2 = sapply(light, run.light, model, 
                     light.scenario = param$light.scenario, diel_group_fw = param$diel_group_fw)
  res2 = cbind.data.frame(t(exts.t2), light, param$t2, param$light.scenario, param$S, param$rep)
  names(res2) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x",
                  "light", "temperature", "light.scenario", "S", "replicate")
  sink(file = 'aaaa',append = T)
  print(names(res1))
  print(dim(res1))
  print(names(res2))
  print(dim(res2))
  sink()
  res = rbind.data.frame(res1,res2)
  # fname = paste("share/groups/tib/ALAN_food_webs/results_ind/", as.character(param$id), sep = '')
  # write.table(res, fname, row.names = F, col.names = F, sep = ",")
  # write.csv(res, fname, row.names = F, col.names = F)
  return(res)
}


############ scripts #################



########  sequential approach ###################
param <- params[[62]]

test <- run.light.gradient(param)

ggplot(test, aes(x= light, y = tot_ext, color = as.factor(temperature))) +
  geom_smooth(method = "lm")


xx = replicate(10, run.light.gradient(masses, n_basal,
                                      t1 = t1, t2 =  t2, S = S,
                                      light.scenario = 'N'),
               simplify = FALSE)

yy = do.call(rbind, xx)
names(yy) = c('exts', 'light', 'temperature')


#########  generate parameter list: ###################

# general params related to food web
reps = 100
temps = c(15, 20)
S.all = c(1, 2)
n_species <- 60
n_basal <- 20
n_cons = n_species - n_basal
n_nut <- 2

# ALAN params
scenarios = c('N', 'C', 'C.N', 'N.D') # light pollution scenarios, which groups are affected
diel_groups = c('D', 'N', 'C') # possible diel groups
probs = c(1,1,1) # probability for a species belonging to each diel group, currently equal probability for all three groups
diel_group_fw = sample(diel_groups, n_cons, replace = TRUE, prob = probs) # associated diel group for each species in the food web

# put everything in a list, 
# each elements contains all parameters for a simulation
params = list()
n = 0
for (i in 1:reps){
  masses <- 10 ^ c(sort(runif(n_basal, 0, 3)),
                   sort(runif(n_species - n_basal, 2, 5)))
  biomasses <- runif(n_species + n_nut, 2, 3)
  diel_group_fw = sample(diel_groups, n_cons, replace = TRUE, prob = probs)
  for (scenario in scenarios){
    for (S in S.all){
      n = n+1
      elem = list(light.scenario = scenario, rep = i, t1 = temps[1], t2 = temps[2], 
                  S = S, masses = masses, biomasses = biomasses,
                  n_species = n_species, n_basal = 20, n_nut = 2,
                  diel_group_fw = diel_group_fw, id = n)
      params[[n]] = elem
    }
  }
}

params[[1]]
################### run them all ################################

# wrapper function as there is no proper replicate function in parallel
wrapper.run.light.gradient = function(param){

  res = run.light.gradient(param)
  # print msg when a replicate is done, to investigate what takes so much time.
  # system(sprintf('echo "%s"', paste0('rep', param$rep, 'done', collapse="", sep = ' ')))
  return(res)
}
# test:
# wrapper.run.light.gradient(params[[2]])
# run:
nb.cores = 25
params2 = params[1:4]
xx = mclapply(params, FUN = wrapper.run.light.gradient,
              mc.cores = nb.cores)
yy = do.call(rbind, xx)

write.csv(yy, '~/share/groups/tib/ALAN_food_webs/results.csv', row.names = F)

# write.table(yy, "~/share/groups/tib/ALAN_food_webs/results.csv",
#             sep = ",", col.names = !file.exists("~/share/groups/tib/ALAN_food_webs/results.csv"), 
#             append = T, row.names = F)

############ plots #################


tab = read.csv("~/share/groups/tib/ALAN_food_webs/results_warm.csv", header = T)
# tabh = read.csv("~/share/groups/tib/ALAN_food_webs/results_high_temps.csv", header = T)
# tab = rbind.data.frame(tab, tabh)
tab$temperature = as.factor(tab$temperature)

# 
# ggplot(tab, aes(x = light, y = ext_night, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.scenario)
#   facet_grid(vars(light.scenario), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_scenario_night.pdf", width = 30, height = 30, units = "cm")
# 
# 
# ggplot(tab, aes(x = light, y = ext_cresp, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.scenario)
#   facet_grid(vars(light.scenario), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_scenario_cresp.pdf", width = 30, height = 30, units = "cm")
# 
# ggplot(tab, aes(x = light, y = ext_day, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.scenario)
#   facet_grid(vars(light.scenario), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_scenario_day.pdf", width = 30, height = 30, units = "cm")
# 
# 
# ggplot(tab, aes(x = light, y = ext_basals, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.scenario)
#   facet_grid(vars(light.scenario), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_scenario_basals.pdf", width = 30, height = 30, units = "cm")
# 
# 
# 
# ggplot(tab, aes(x = light, y = tot_ext, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.scenario)
#   facet_grid(vars(light.scenario), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_scenario_total.pdf", width = 30, height = 30, units = "cm")
# 
# 

library(reshape2)
# remove total extinctions (problem of scale in facet wraps)
tabx = tab[,-1]

# for S = 1
tab1 = melt(subset(tabx, S ==1), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.scenario', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)

ggplot(tab1, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.scenario)
  facet_grid(vars(light.scenario), vars(Species))+
  ggtitle(paste("S= ", 1, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S1.pdf", width = 30, height = 30, units = "cm")


tab2 = melt(subset(tabx, S ==2), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.scenario', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab2, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.scenario)
  facet_grid(vars(light.scenario), vars(Species))+
  ggtitle(paste("S= ", 2, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S2.pdf", width = 30, height = 30, units = "cm")

tab5 = melt(subset(tabx, S ==5), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.scenario', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab5, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.scenario)
  facet_grid(vars(light.scenario), vars(Species))+
  ggtitle(paste("S= ", 5, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S5.pdf", width = 30, height = 30, units = "cm")


tab20 = melt(subset(tabx, S ==20), 
             measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
             id.vars = c('x', 'light', 'temperature', 'light.scenario', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab20, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.scenario)
  facet_grid(vars(light.scenario), vars(Species))+
  ggtitle(paste("S= ", 20, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S20.pdf", width = 30, height = 30, units = "cm")



