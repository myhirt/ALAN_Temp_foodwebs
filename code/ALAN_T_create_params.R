library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

set.seed(123)

#########  generate parameter list: ###################
fws = 100
temps = c(15, 20)
S.all = c(1)
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
for (i in 1:fws){
  masses <- 10 ^ c(sort(runif(n_basal, 0, 3)),
                   sort(runif(n_species - n_basal, 2, 5)))
  biomasses <- runif(n_species + n_nut, 2, 3)
  period = sample(types, n_cons, replace = TRUE, prob = probs)
    for (S in S.all){
      n = n+1
      elem = list(fw = i, t1 = temps[1], t2 = temps[2], 
                  S = S, masses = masses, biomasses = biomasses,
                  n_species = n_species, n_basal = 20, n_nut = 2,
                  period = period)
      params[[n]] = elem
  }
}

# Create a function to extract and flatten each element
flatten_param <- function(p) {
  tibble(
    fw = p$fw,
    t1 = p$t1,
    t2 = p$t2,
    S = p$S,
    n_species = p$n_species,
    n_basal = p$n_basal,
    n_nut = p$n_nut,
    masses = paste(round(p$masses, 4), collapse = ";"),
    biomasses = paste(round(p$biomasses, 4), collapse = ";"),
    period = paste(p$period, collapse = ";")
  )
}

# Apply to all elements and bind into one data frame
params_df <- map_dfr(params, flatten_param)
as.numeric(strsplit(params_df$masses[1], ";")[[1]])

# Save to CSV
write.csv(params_df, "HPC_ALAN_T_params.csv", row.names = FALSE)
