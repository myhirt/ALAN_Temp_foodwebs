library(ATNr, lib.loc='~/share/groups/tib/ALAN_food_webs/')
library(ggplot2)
library(parallel)
set.seed(12)

# force BLAS to use only one core
library('RhpcBLASctl')
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

rm(list = ls())

############ functions #################

overlap.modification = function(x, act.time, effect = "control"){
  # update the "time overlap matrix"
  # activity time is the vector containing consumers' activity period
  # effect controls which type of species light has an effect on:
    # "control": no effect
    # "N": light only affects nocturnal species
    # "C": light only affects crepuscular species
    # "CN": light affects both nocturnal and crepuscular species
  
  mat.act = outer(act.time, act.time, paste, sep = ".")
  mat.overlap = matrix(0, nrow = length(act.time), ncol = length(act.time))
  
  # control
  mat.overlap[mat.act == "D.D" | mat.act == "N.N"] = 1
  mat.overlap[grep('C', mat.act)] = 0.5
  
  # light only affects nocturnal species
  if (effect == "N"){
    mat.overlap[mat.act == "N.N"] = 1
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 - x
  }
  # only affects crepuscular species
  if (effect == "C"){
    mat.overlap[mat.act == "C.D" | mat.act == "D.C"] = 0.5 - x
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 + x
  }
  # affects crepuscular and nocturnal species
  if (effect == "C.N" | effect == "N.C"){
    mat.overlap[mat.act == "C.D" | mat.act == "D.C"] = 0.5 - x
    # mat.overlap[mat.act == "C.N" | mat.act == "N.C" | mat.act == "C.C"] = 0.5 + x/2
    # mat.overlap[mat.act == "N.N"] = 1 - x/2
    
  } 
  # test scenarion: night species shift to day
  if (effect == "N.D"){
    mat.overlap[mat.act == "D.N" | mat.act == "N.D"] = 0.0 + x
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 + x
    # mat.overlap[mat.act == "N.N"] = 1 - x/2
    
  } 
  
  # mat.1s = matrix(1, nrow = n_basal, ncol = n_species - n_basal)
  # mat.overlap = rbind(mat.1s, mat.overlap)
  return(mat.overlap)
}

run.light = function(x, model, light.effect, period){
  # This function returns the number of extinctions of a given model 
  # associated with a light intensity x
  # light effect is the scenario: which species types are affected by light pollution
  # period is a vector containing the activity period of non basal species
  
  # modifications due to activity overlap
  overlaps = overlap.modification(x, period, light.effect)
  # update b accordingly
  b1 = model$b
  model$b[(n_basal +1):n_species, ] = model$b[(n_basal +1):n_species, ] * overlaps
  model$w[model$w!=0] = 1
  new.basals = sum(colSums(model$b) == 0)
  if (new.basals > 0){print(new.basals)}
  # running dynamics
  times <- seq(0, 100000, 100)
  
  sol = tryCatch(
    {
      R.utils::withTimeout({
    
        time.series = lsoda_wrapper(times, biomasses, model, verbose = FALSE)
        extinctions <- time.series[nrow(time.series), (n_nut+2):ncol(time.series)] < 1e-6
        bioms = time.series[nrow(time.series), (n_nut+2):ncol(time.series)]
        animal.bioms = tail(bioms, n_cons)
        basal.ext = sum(extinctions[1:n_basal])
        basal.ext = (model$nb_b - basal.ext)/model$nb_b 
        tot.ext = sum(!extinctions)/model$nb_s
        animal.ext = tail(extinctions, n_cons)
        night.ext = sum(!animal.ext[period == 'N'])/ sum(period == 'N')
        cresp.ext = sum(!animal.ext[period == 'C'])/ sum(period == 'C')
        day.ext = sum(!animal.ext[period == 'D'])/ sum(period == 'D')
        basal.bioms = sum(bioms[1:n_basal])
        night.biom = sum(animal.bioms[period == "N"])
        cresp.biom = sum(animal.bioms[period == 'C'])
        day.biom = sum(animal.bioms[period == 'D']) 
        model$b = b1
        to.return = c(tot.ext, basal.ext, night.ext, cresp.ext, day.ext, 
                      basal.bioms, night.biom, cresp.biom, day.biom, x)
        
        names(to.return) = c("tot.ext", "ext.basals", "ext.night", "ext.cresp", "ext.day", 
                             "basal_bioms", "night_biom", "cresp_biom", "day_.biom", "x")
        return(to.return)
      }, 
      timeout = 10)
    }, 
    TimeoutException  = function(x){
      return(rep(NA, 10))
    },
    warning = function(w){return(rep(NA, 10))}
  )
  
  
#   sol <- R.utils::withTimeout(
#     ,
#     timeout = 8,
#     onTimeout = "silent"
#   )
#   extinctions <- sum(sol[nrow(sol), (n_nut+2):ncol(sol)] < 1e-6)
#   return(extinctions)
}

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
  # even without light effect
  # if yes, cancel run
  overlaps = overlap.modification(x, period, 0)
  b2 = model$b
  b2[(n_basal +1):n_species, ] = b2[(n_basal +1):n_species, ] * overlaps
  if (any(colSums(b2) == 0.0)){return(NULL)}
  
  exts.t1 = sapply(light, run.light, model, 
                     light.effect = param$light.effect, period = param$period)
  res1 = cbind.data.frame(t(exts.t1), light, param$t1, param$light.effect, param$S, param$rep)
  names(res1) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x", 
                  "light", "temperature", "light.effect", "S", "replicate")
  
  model <- initialise_default_Unscaled_nuts(model, L, temperature = param$t2)
  model$S = rep(param$S, param$n_nut)
  model$q = rep(1.2, n_species - n_basal)
  exts.t2 = sapply(light, run.light, model, 
                     light.effect = param$light.effect, period = param$period)
  res2 = cbind.data.frame(t(exts.t2), light, param$t2, param$light.effect, param$S, param$rep)
  names(res2) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x",
                  "light", "temperature", "light.effect", "S", "replicate")
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
# xx = replicate(10, run.light.gradient(masses, n_basal, 
#                                       t1 = t1, t2 =  t2, S = S, 
#                                       light.effect = 'N'), 
#                simplify = FALSE)
# yy = do.call(rbind, xx)
# names(yy) = c('exts', 'light', 'temperature')


#########  generate parameter list: ###################


reps = 100
effects = c('N', 'C', 'C.N', 'N.D')
temps = c(15, 20)
S.all = c(1, 2)
n_species <- 60
n_basal <- 20
n_cons = n_species - n_basal
n_nut <- 2

types = c('D', 'N', 'C')
probs = c(1,1,1) # if we want more from certain types
period = sample(types, n_cons, replace = TRUE, prob = probs)

# put everything in a list, 
# each elements contains all parameters for a simulation
params = list()
n = 0
for (i in 1:reps){
  masses <- 10 ^ c(sort(runif(n_basal, 0, 3)),
                   sort(runif(n_species - n_basal, 2, 5)))
  biomasses <- runif(n_species + n_nut, 2, 3)
  period = sample(types, n_cons, replace = TRUE, prob = probs)
  for (effect in effects){
    for (S in S.all){
      n = n+1
      elem = list(light.effect = effect, rep = i, t1 = temps[1], t2 = temps[2], 
                  S = S, masses = masses, biomasses = biomasses,
                  n_species = n_species, n_basal = 20, n_nut = 2,
                  period = period, id = n)
      params[[n]] = elem
    }
  }
}

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
#   # facet_wrap(~ light.effect)
#   facet_grid(vars(light.effect), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_effect_night.pdf", width = 30, height = 30, units = "cm")
# 
# 
# ggplot(tab, aes(x = light, y = ext_cresp, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.effect)
#   facet_grid(vars(light.effect), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_effect_cresp.pdf", width = 30, height = 30, units = "cm")
# 
# ggplot(tab, aes(x = light, y = ext_day, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.effect)
#   facet_grid(vars(light.effect), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_effect_day.pdf", width = 30, height = 30, units = "cm")
# 
# 
# ggplot(tab, aes(x = light, y = ext_basals, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.effect)
#   facet_grid(vars(light.effect), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_effect_basals.pdf", width = 30, height = 30, units = "cm")
# 
# 
# 
# ggplot(tab, aes(x = light, y = tot_ext, color = temperature))+
#   # geom_point(alpha = 0.1, shape = 20)+
#   stat_smooth()+
#   # facet_wrap(~ light.effect)
#   facet_grid(vars(light.effect), vars(S))
# ggsave("~/share/groups/tib/ALAN_food_webs/light_effect_total.pdf", width = 30, height = 30, units = "cm")
# 
# 

library(reshape2)
# remove total extinctions (problem of scale in facet wraps)
tabx = tab[,-1]

# for S = 1
tab1 = melt(subset(tabx, S ==1), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.effect', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)

ggplot(tab1, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.effect)
  facet_grid(vars(light.effect), vars(Species))+
  ggtitle(paste("S= ", 1, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S1.pdf", width = 30, height = 30, units = "cm")


tab2 = melt(subset(tabx, S ==2), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.effect', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab2, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.effect)
  facet_grid(vars(light.effect), vars(Species))+
  ggtitle(paste("S= ", 2, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S2.pdf", width = 30, height = 30, units = "cm")

tab5 = melt(subset(tabx, S ==5), 
            measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
            id.vars = c('x', 'light', 'temperature', 'light.effect', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab5, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.effect)
  facet_grid(vars(light.effect), vars(Species))+
  ggtitle(paste("S= ", 5, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S5.pdf", width = 30, height = 30, units = "cm")


tab20 = melt(subset(tabx, S ==20), 
             measure.vars = c('pers_basals', 'pers_night', 'pers_cresp', 'pers_day'),
             id.vars = c('x', 'light', 'temperature', 'light.effect', 'S', 'replicate'),
            variable.name = 'Species',
            value.name = 'Persistence'
)
ggplot(tab20, aes(x = light, y = Persistence, color = temperature))+
  # geom_point(alpha = 0.1, shape = 20)+
  stat_smooth()+
  # facet_wrap(~ light.effect)
  facet_grid(vars(light.effect), vars(Species))+
  ggtitle(paste("S= ", 20, sep = ''))
ggsave("~/share/groups/tib/ALAN_food_webs/S20.pdf", width = 30, height = 30, units = "cm")



