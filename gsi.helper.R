##### Function for assisting in calculating the genealogical sortng index (gsi) #####

gsi.helper <- function(spal_matrix, population_set, gene_trees){ # spal_matrix = species-allele matrix, population_set = names of populations, gene_trees = trees to be analyzed
  
  gsi_locus <- c() # create vector for storing gsi values for all trees
  
  for (g in gene_trees){ # for each predictive tree
    
    gsi_vector <- c() # create vector for storing gsi values for each pop
    
    for (pop in population_set){ # for each population
      grpin <- which(spal_matrix[,2] == pop) # get indices for alleles belonging to population
      grp <- g$tip.label[grpin] # get tip labels for alleles belonging to population
      gsi_vector <- c(gsi_vector, gsi.calc(g, grp)) # calculate gsi for population and add to results list
    }
    
    gsi_mean <- mean(as.numeric(gsi_vector)) # average gsi values
    gsi_locus <- c(gsi_locus, gsi_mean) # add gsi value to locus results
    
  }
  
  return(gsi_locus)
  
}



