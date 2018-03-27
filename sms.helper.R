##### Function for assisting in calculating Slatkin and Maddison's s #####

sms.helper <- function(spal_matrix, gene_trees){ # spal_matrix = species-allele matrix, population_set = names of populations, gene_trees = trees to be analyzed
  
  sms_locus <- c() # create vector for storing sms values for all trees
  
  for (g in gene_trees){ # for each predictive tree
    
    sms_locus <- c(sms_locus, sms.calc(g, spal_matrix)) # calculate sms and add value to locus results
    
  }
  
  return(sms_locus)
  
}