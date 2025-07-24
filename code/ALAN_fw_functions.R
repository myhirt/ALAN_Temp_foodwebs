# This function constructs a time-overlap matrix for species in a food web, based on their diel activity patterns and a 
# specified light pollution scenario. This matrix can then be used to modulate feeding interactions 
# according to how much overlap there is in the times species are active.

overlap.modification = function(x, diel.groups, light.scenario = "control"){
  # "x" is a parameter (between 0 and 0.5) that determines how much the overlap is modified by light pollution
  # " diel.groups" is the vector containing consumers' diel group
  # "light.scenario" refers to which species group is affected by light pollution
  #     "control": no effect
  #     "N": light only affects nocturnal species
  #     "C": light only affects crepuscular species
  #     "CN": light affects both nocturnal and crepuscular species
  #     "CD": light affects both crepuscular and diurnal species
  
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
  if (light.scenario == "C.D"){
    mat.overlap[mat.diel == "D.N" | mat.diel == "N.D"] = 0.0 + x
    mat.overlap[mat.diel == "C.N" | mat.diel == "N.C"] = 0.5 + x
    # mat.overlap[mat.diel == "C.N" | mat.diel == "N.C" | mat.diel == "C.C"] = 0.5 + x/2
    # mat.overlap[mat.diel == "N.N"] = 1 - x/2
    
  } 
  
  # mat.1s = matrix(1, nrow = n_basal, ncol = n_species - n_basal)
  # mat.overlap = rbind(mat.1s, mat.overlap)
  return(mat.overlap)
}


#The sim.ALAN function is a wrapper for running food web simulations under different light pollution 
# intensities (x = 0-0.5) and diel niche disturbance scenarios (light.scenario)
# 
# It outputs
#   Total extinction
#   Group-specific extinctions (night, crepuscular, day)
#   Biomass remaining in each group


sim.ALAN <- function(x, masses, biomasses, n_basal, n_nut,  probs, light.scenario, diel_group, temp) {
  # x: overlap modification intensity (0 to 0.5)
  # light.scenario: diel groups affected by light pollution
  # diel_group_fw: diel group (N, D, C) of each consumer species
  # Create L matrix and food web structure
  L <- create_Lmatrix(masses, n_basal, Ropt = 100, gamma = 2, th = 0.01)
  fw <- (L > 0) * 1  # binary adjacency matrix
  
  model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
  model <- initialise_default_Unscaled_nuts(model, L, temperature = temp)
  model$S = rep(2, n_nut)
  
  # Modify consumer-consumer interaction strengths based on diel overlap
  overlaps <- overlap.modification(x, diel_group, light.scenario)
  model$b[(n_basal + 1):n_species, ] <- model$b[(n_basal + 1):n_species, ] * overlaps
  model$w[model$w != 0] <- 1
  new.basals = sum(colSums(model$b) == 0)
  if (new.basals > 0){print(new.basals)}
  
  # Exit early if modification disconnects consumers from all basal resources
  if (any(colSums(model$b) == 0)) return(NULL)
  
  times <- seq(0, 1000, 100)
  
  # Run dynamics with timeout protection
  tryCatch({
    R.utils::withTimeout({
      ts <- lsoda_wrapper(times, biomasses, model, verbose = FALSE)
      final <- ts[nrow(ts), (n_nut + 2):ncol(ts)]
      extinct <- final < 1e-6
      animal_biom <- tail(final, n_cons)
      
      out <- data.frame(
        tot.pers   = sum(!extinct) / model$nb_s,
        basal.pers = (model$nb_b - sum(extinct[1:n_basal])) / model$nb_b,
        night.pers = mean(!tail(extinct, n_cons)[diel_group == 'N']),
        crep.pers = mean(!tail(extinct, n_cons)[diel_group == 'C']),
        day.pers   = mean(!tail(extinct, n_cons)[diel_group == 'D']),
        basal.bioms = sum(final[1:n_basal]),
        night.biom  = sum(animal_biom[diel_group == 'N']),
        cresp.biom  = sum(animal_biom[diel_group == 'C']),
        day.biom    = sum(animal_biom[diel_group == 'D']),
        x = x,
        light.scenario = light.scenario,
        temp = temp
      )
      return(out)
    }, timeout = 100)
  },
  TimeoutException = function(e) rep(NA, 12),
  warning = function(w) rep(NA, 12))
}







