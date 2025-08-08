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

run.light.gradient = function(light.effect, t1, t2, S, n_basal, n_species, n_nuts, masses, biomasses, period, rep){
  # run a light intensity gradient at 2 temperatures, 
  n_cons = n_species - n_basal
  # create the L matrix
  L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)
  # create the 0/1 version of the food web
  fw <- L
  fw[fw > 0] <- 1
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  
  light = seq(0, 0.5, by = 0.05)
  model <- initialise_default_Unscaled_nuts(model, L, temperature = t1)
  model$S = rep(S, n_nut)
  model$q = rep(1.2, n_species - n_basal)
  
  # test if the separation btween nocturnal and diurnal lead to some new basal species
  # even without light effect
  # if yes, cancel run
  overlaps = overlap.modification(x, period, 0)
  b2 = model$b
  b2[(n_basal +1):n_species, ] = b2[(n_basal +1):n_species, ] * overlaps
  if (any(colSums(b2) == 0.0)){return(NULL)}
  
  exts.t1 = sapply(light, run.light, model, 
                   light.effect = light.effect, period = period)
  res1 = cbind.data.frame(t(exts.t1), light, t1, light.effect, S, fw)
  names(res1) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x", 
                  "light", "temperature", "light.effect", "S", "fw")
  
  model <- initialise_default_Unscaled_nuts(model, L, temperature = t2)
  model$S = rep(S, n_nut)
  model$q = rep(1.2, n_species - n_basal)
  exts.t2 = sapply(light, run.light, model, 
                   light.effect = light.effect, period = period)
  res2 = cbind.data.frame(t(exts.t2), light, t2, light.effect, S, fw)
  names(res2) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x",
                  "light", "temperature", "light.effect", "S", "fw")
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
