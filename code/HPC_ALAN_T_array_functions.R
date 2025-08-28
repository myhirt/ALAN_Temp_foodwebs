
overlap.modification = function(x, act.time, effect = "control"){
  # update the "time overlap matrix"
  mat.act = outer(act.time, act.time, paste, sep = ".")
  mat.overlap = matrix(0, nrow = length(act.time), ncol = length(act.time))
  
  # baseline (control, no light effect)
  mat.overlap[mat.act == "D.D" | mat.act == "N.N"] = 1
  mat.overlap[grep('C', mat.act)] = 0.5
  
  # nocturnal species affected
  if (effect == "N"){
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 - x
  }
  # crepuscular species affected
  if (effect == "C"){
    mat.overlap[mat.act == "C.D" | mat.act == "D.C"] = 0.5 - x
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 + x
  }
  # nocturnal + crepuscular affected
  if (effect %in% c("C.N","N.C")){
    mat.overlap[mat.act == "C.D" | mat.act == "D.C"] = 0.5 - x
  }
  # nocturnal shift to diurnal
  if (effect == "N.D"){
    mat.overlap[mat.act == "D.N" | mat.act == "N.D"] = 0.0 + x
    mat.overlap[mat.act == "C.N" | mat.act == "N.C"] = 0.5 + x
  }
  
  return(mat.overlap)
}


run.light = function(x, model, light.effect, period){
  # This function returns the number of extinctions of a given model 
  # associated with a light intensity x
  # light effect is the scenario: which species types are affected by light pollution
  # period is a vector containing the activity period of non basal species
  
  # modifications due to activity overlap
  overlaps = overlap.modification(x, period, light.effect)
  
  # update b matrix
  b1 = model$b
  model$b[(n_basal +1):n_species, ] = model$b[(n_basal +1):n_species, ] * overlaps
  model$w[model$w!=0] = 1
  
  times <- seq(0, 100000, 100)   # unified time horizon
  
  sol = tryCatch(
    {
      R.utils::withTimeout({
        time.series = lsoda_wrapper(times, biomasses, model, verbose = FALSE)
        extinctions <- time.series[nrow(time.series), (n_nut+2):ncol(time.series)] < 1e-6
        bioms = time.series[nrow(time.series), (n_nut+2):ncol(time.series)]
        animal.bioms = tail(bioms, n_cons)
        
        # persistence fractions
        basal.pers = (model$nb_b - sum(extinctions[1:n_basal]))/model$nb_b 
        tot.pers   = sum(!extinctions)/model$nb_s
        animal.ext = tail(extinctions, n_cons)
        night.pers = sum(!animal.ext[period == 'N'])/ sum(period == 'N')
        cresp.pers = sum(!animal.ext[period == 'C'])/ sum(period == 'C')
        day.pers   = sum(!animal.ext[period == 'D'])/ sum(period == 'D')
        
        # biomasses
        basal.bioms = sum(bioms[1:n_basal])
        night.biom  = sum(animal.bioms[period == "N"])
        cresp.biom  = sum(animal.bioms[period == 'C'])
        day.biom    = sum(animal.bioms[period == 'D']) 
        
        model$b = b1  # reset
        
        to.return = c(tot.pers, basal.pers, night.pers, cresp.pers, day.pers,
                      basal.bioms, night.biom, cresp.biom, day.biom, x)
        names(to.return) = c("pers_tot", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                             "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x")
        return(to.return)
      }, timeout = 10)
    },
    TimeoutException  = function(x){ return(rep(NA, 10)) },
    warning = function(w){ return(rep(NA, 10)) }
  )
}


run.light.gradient = function(light.effect, t1, t2, S, n_basal, n_species, n_nut, masses, biomasses, period, rep){
  # run a light intensity gradient at 2 temperatures, 
  
  # create the L matrix
  L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)
  # create the 0/1 version of the food web
  fw <- L
  fw[fw > 0] <- 1
  connectance = sum(fw)/(n_species*n_species)
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  
  light = seq(0, 0.5, by = 0.05)
  model <- initialise_default_Unscaled_nuts(model, L, temperature = t1)
  model$S = rep(S, n_nut)
  model$q = rep(1.2, n_species - n_basal)
  
  # test if the separation btween nocturnal and diurnal lead to some new basal species
  # even without light effect
  # if yes, cancel run
  overlaps = overlap.modification(x, period, "control")
  b2 = model$b
  b2[(n_basal +1):n_species, ] = b2[(n_basal +1):n_species, ] * overlaps
  if (any(colSums(b2) == 0.0)){return(NULL)}
  
  exts.t1 = sapply(light, run.light, model, 
                   light.effect = light.effect, period = period)
  res1 = cbind.data.frame(t(exts.t1), light, t1, light.effect, S, rep, connectance)
  names(res1) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x", 
                  "light", "temperature", "light.effect", "S", "replicate", "connectance")
  
  model <- initialise_default_Unscaled_nuts(model, L, temperature = t2)
  model$S = rep(S, n_nut)
  model$q = rep(1.2, n_species - n_basal)
  exts.t2 = sapply(light, run.light, model, 
                   light.effect = light.effect, period = period)
  res2 = cbind.data.frame(t(exts.t2), light, t2, light.effect, S, rep, connectance)
  names(res2) = c("tot_ext", "pers_basals", "pers_night", "pers_cresp", "pers_day", 
                  "basal_bioms", "night_biom", "cresp_biom", "day_biom", "x",
                  "light", "temperature", "light.effect", "S", "replicate", "connectance")
  # sink(file = 'aaaa',append = T)
  # print(names(res1))
  # print(dim(res1))
  # print(names(res2))
  # print(dim(res2))
  # sink()
  res = rbind.data.frame(res1,res2)
  # fname = paste("share/groups/tib/ALAN_food_webs/results_ind/", as.character(param$id), sep = '')
  # write.table(res, fname, row.names = F, col.names = F, sep = ",")
  # write.csv(res, fname, row.names = F, col.names = F)
  return(res)
}
