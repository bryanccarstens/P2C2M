#############################
#   CALCULATE SUM STATS   ###
#############################

sum_stats <- function(bayesallfile, mig_tree, loci, samples, pops, post_samples, descrStats, pop_assign, singleAllele, empirical_popgenome){

########## Read in post pred infinite allele data ##########
  cat("\nReading predictive data\n")
  ppt_infsites <- list()
  ppt_infsites <- lapply(seq(1, length(loci), 1), function(x) PopGenome::readMS(sprintf("Locus%s_pps.trees", x)))
  
# Define populations
  pred_pop_defs <- lapply(seq(1, length(loci), 1), function(x) samples_per_locus(samples[,x], pops)) # create list of vectors with each vector containing the individuals in each pop
  ppt_infsites <- mapply(PopGenome::set.populations, ppt_infsites, pred_pop_defs) # define the pops for each post pred dataset
  cat("\n\n")
  
########## Read in migrate gene trees ##########

  trees <- migtree.reader(bayesallfile, mig_tree, post_samples, length(loci)) # use migtree.reader to read in tree subset matching predictive trees
  trees <- lapply(trees, function(x) as.character(x))
  migrate_trees <- ape::read.tree(text = unlist(trees)) # read all gene trees into ape
  cat("\n")
  
  if ("GSI" %in% descrStats | "SMS" %in% descrStats){ # if user selects GSI or SMS statistics
    pred_trees <- list() # create empty list for storing pred trees
    post_trees <- list() # create list for storing post trees
    for (loc in 1:length(loci)){ # for each locus
      print(paste0("Reading predictive trees from Locus", loc))
      file=paste("Locus", loc, "_pps.trees", sep="") # create file name for post pred simulation tree file
      pred_trees[loc] <- list(ape::read.tree(file)) # read in post pred trees
      
      loc_ind <- seq((length(post_samples) * loc) - length(post_samples) + 1, (length(post_samples) * loc), 1) # get indices of trees in tree list for the current locus
      post_trees[loc] <- list(migrate_trees[loc_ind]) # subset trees for the current locus
      
    }
  }
  cat("\n")
  
########## Calculate summary statistics ##########

  spal <- spal.creator(pop_assign, singleAllele) # manipulate species-allele matrix and create matrix for posterior predictive trees
  
  stats <- calc.sum.stats(ppt_infsites, empirical_popgenome, pred_trees, post_trees, post_samples, spal, loci, descrStats)
 
  
  return(stats) 
}

