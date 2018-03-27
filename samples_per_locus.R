##### Function to produce list of alleles in each population for posterior predictive simulation datasets #####

samples_per_locus <- function(samples_per_loc, pops){

  inds_per_pop <- c(0)
  inds_per_loc <- list()
  for (p in 1:pops){
    inds_per_pop <- c(inds_per_pop, sum(samples_per_loc[c(1:p)])) # calculate cumulative number of individuals for population
    inds <- seq(inds_per_pop[p] + 1, inds_per_pop[p + 1], 1) # create individual numbers for population
    inds_per_loc[p] <- list(inds)
  }
  
  return(inds_per_loc)
}
