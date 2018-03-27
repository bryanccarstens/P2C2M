##### Function for modifying the species allele matrix and identifying population sets for the summary stats #####

spal.creator <- function(pop_assign, singleAllele){ # pop_assign = current species-allele matrix; singleAllele = any species (populations) containing only a single individual
  
  ### Create species allele matrices for post and post pred ###
  post_spal <- do.call("rbind", pop_assign) # create a combined species (population)-allele matrix
  pred_spal <- post_spal
  al_list <- lapply(unique(pred_spal[,2]), function(x) which(pred_spal[,2] == x)) # create lists of species allele associations
  pred_spal[,1] <- unlist(al_list)
  
  ### Create population sets for post and post pred ###
  post_popset <- unique(post_spal[,2]) # get unique population names
  pred_popset <- unique(pred_spal[,2]) # get unique population numbers
  
  if (length(singleAllele) > 0){ # if any population consists of only a single individual
    post_popset_nosingles <- suppressWarnings(post_popset[post_popset != singleAllele]) # remove single allele populations (for gsi)
    pred_singleAllele <- which(post_popset == singleAllele) # get index of single allele population
    pred_popset_nosingles <- suppressWarnings(pred_popset[pred_popset != pred_singleAllele]) # remove single allele populations
  }
  
  else {
    post_popset_nosingles <- post_popset
    pred_popset_nosingles <- pred_popset
  }
  
  spal <- list()
  spal$post <- post_spal
  spal$pred <- pred_spal
  spal$post_popset <- post_popset
  spal$pred_popset <- pred_popset
  spal$post_popset_nosingles <- post_popset_nosingles
  spal$pred_popset_nosingles <- pred_popset_nosingles
  
  return(spal)
}
