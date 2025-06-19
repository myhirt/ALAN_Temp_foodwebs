# This function constructs a time-overlap matrix for species in a food web, based on their diel activity patterns and a 
# specified light pollution light.scenario. This matrix can then be used to modulate feeding interactions 
# according to how much overlap there is in the times species are active.

overlap.modification = function(x, diel.groups, light.scenario = "control"){
  # "x" a parameter (between 0 and 0.5) that determines how much the overlap is modified by light pollution
  # " diel.groups" is the vector containing consumers' diel group
  # "light.scenario" refers to the light.scenario which species group is affected by light pollution
  # "control": no effect
  # "N": light only affects nocturnal species
  # "C": light only affects crepuscular species
  # "CN": light affects both nocturnal and crepuscular species
  
  # matrix that contains all possible combinations of interacting diel groups
  mat.diel = outer(diel.groups, diel.groups, paste, sep = ".")
  # creates empty matrix that will be filled with the overlap values between the diel groups based on the light.scenario
  mat.overlap = matrix(0, nrow = length(diel.groups), ncol = length(diel.groups))
  
  # control ("normal overlap", e.g. diurnal interatc with diurnal and crep. but not nocturnal)
  mat.overlap[mat.diel == "D.D" | mat.diel == "N.N"] = 1
  mat.overlap[grep('C', mat.diel)] = 0.5
  
  # light only affects nocturnal species
  if (light.scenario == "N"){
    mat.overlap[mat.diel == "C.N" | mat.diel == "N.C"] = 0.5 - x # overlap between nocturnal and crepuscular species is decreased by value x
  }
  # only affects crepuscular species
  if (light.scenario == "C"){
    mat.overlap[mat.diel == "C.D" | mat.diel == "D.C"] = 0.5 - x
    mat.overlap[mat.diel == "C.N" | mat.diel == "N.C"] = 0.5 + x
  }
  # affects crepuscular and nocturnal species
  if (light.scenario == "C.N" | light.scenario == "N.C"){
    mat.overlap[mat.diel == "C.D" | mat.diel == "D.C"] = 0.5 - x
    # mat.overlap[mat.diel == "C.N" | mat.diel == "N.C" | mat.diel == "C.C"] = 0.5 + x/2
    # mat.overlap[mat.diel == "N.N"] = 1 - x/2
    
  } 
  # test light.scenarion: night species shift to day
  if (light.scenario == "N.D"){
    mat.overlap[mat.diel == "D.N" | mat.diel == "N.D"] = 0.0 + x
    mat.overlap[mat.diel == "C.N" | mat.diel == "N.C"] = 0.5 + x
    # mat.overlap[mat.diel == "N.N"] = 1 - x/2
    
  } 
  
  # mat.1s = matrix(1, nrow = n_basal, ncol = n_species - n_basal)
  # mat.overlap = rbind(mat.1s, mat.overlap)
  return(mat.overlap)
}


#The sim.ALAN function is a wrapper for running food web simulations under different light pollution 
# intensities (x = 0-0.5) and diel niche disturbance light.scenarios (light.scenario)
# 
# It outputs key results like:
#   Total extinction
#   Group-specific extinctions (night, crespuscular, day)
#   Biomass remaining in each group


sim.ALAN = function(x, masses, n_basal, light.scenario, diel_group_fw, temp) {
  # "x" a parameter (between 0 and 0.5) that determines how much the overlap is modified by light pollution
  # light.scenario: which species diel.groups are affected by light pollution
  # diel_group_fw is a vector containing the diel group of consumer species in the food web 
  
  # create the L matrix
  L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)
  # create the 0/1 version of the food web
  fw <- L
  fw[fw > 0] <- 1
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  model <- initialise_default_Unscaled_nuts(model, L, temperature = temp)
  # calculate overlap in activity time due to diel group and light pollution intensity (x)
  overlaps = overlap.modification(x, diel_group_fw, light.scenario)
  # update attack rate b between consumer species accordingly
  b1 = model$b 
  model$b[(n_basal +1):n_species, ] = model$b[(n_basal +1):n_species, ] * overlaps # modify based on overlap
  model$w[model$w!=0] = 1
  #test if the separation btween nocturnal and diurnal lead to some new basal species
  if (any(colSums(model$b) == 0.0)){return(NULL)}
  
  # running dynamics
  times <- seq(0, 1000, 100)
  
  sol = tryCatch(
    {
      R.utils::withTimeout({
        
        time.series = lsoda_wrapper(times, biomasses, model, verbose = FALSE)
        extinctions <- time.series[nrow(time.series), (n_nut+2):ncol(time.series)] < 1e-6
        bioms = time.series[nrow(time.series), (n_nut+2):ncol(time.series)]
        animal.bioms = tail(bioms, n_cons)
        basal.ext = sum(extinctions[1:n_basal])
        basal.ext = (model$nb_b - basal.ext)/model$nb_b 
        tot.ext = sum(extinctions)/model$nb_s
        animal.ext = tail(extinctions, n_cons)
        night.ext = sum(!animal.ext[diel_group_fw == 'N'])/ sum(diel_group_fw == 'N')
        cresp.ext = sum(!animal.ext[diel_group_fw == 'C'])/ sum(diel_group_fw == 'C')
        day.ext = sum(!animal.ext[diel_group_fw == 'D'])/ sum(diel_group_fw == 'D')
        basal.bioms = sum(bioms[1:n_basal])
        night.biom = sum(animal.bioms[diel_group_fw == "N"])
        cresp.biom = sum(animal.bioms[diel_group_fw == 'C'])
        day.biom = sum(animal.bioms[diel_group_fw == 'D']) 
        to.return = as.data.frame(cbind(tot.ext, basal.ext, night.ext, cresp.ext, day.ext, 
                      basal.bioms, night.biom, cresp.biom, day.biom, x))
        to.return$light.scenario <- light.scenario
        to.return$x <- x
        to.return$temp <- temp
        return(to.return)
      }, 
      timeout = 10)
    }, 
    TimeoutException  = function(x){
      return(rep(NA, 10))
    },
    warning = function(w){return(rep(NA, 10))}
  )
}







