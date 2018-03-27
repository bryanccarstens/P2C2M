##### Function to read in and format gene trees from migrate output #####

migtree.reader <- function(bayesallfile, treefile, post_samples, n_locs){ # Function for reading genetrees from migrate treefiles. post_samples = indices of trees that match those subset in the bayesallfile; nlocs = number of loci
  post_samples2 <- post_samples
  post_samples3 <- c(0, post_samples2[1:length(post_samples2)-1]) + 1 # create vector to be used below
  post_samples4 <- post_samples2 - post_samples3 # subtract the two vectors to obtain a vector of values, each being the number of trees to skip
  post_samples5 <- (post_samples4 * 3) # convert the number of trees to skip to the number of lines to skip
  post_samples5[1] <- (post_samples2[1] * 3) # change the first value to account for locus 1 only having 1 more tree than the bayesallfile instead of two
  post_samples5 <- post_samples5 + 2 # correct line numbers to get only the lines with trees
  leftover <- ((nrow(bayesallfile[bayesallfile$Locus==1,]) - post_samples[length(post_samples)]) * 3) - 2 # calculate the remaining number of lines to skip between the last tree of one locus and first tree of the next locus

  file <- file(treefile, "r")
  mig_trees <- c()
  
  print("Reading posterior trees from Locus1")
  for (samp in post_samples5){
    readLines(file, n=samp) # skip lines
    mig_trees <- c(mig_trees, readLines(file, n=1)) # read tree
  }
  
  for (y in 2:(n_locs)){
    print(paste0("Reading posterior trees from Locus", y))
    readLines(file, n=(5 + leftover)) # skip lines to first tree of next locus
    for (samp in post_samples5){
      readLines(file, n=samp) # skip lines
      mig_trees <- c(mig_trees, readLines(file, n=1)) # read tree
    }
  }
  
  close(file)
  
  mig_trees <- gsub("\\s+\\[(.*?)\\]", "", mig_trees) # get rid of migration estimates within the gene trees
  mig_trees <- gsub("=", "", mig_trees) # remove "=" from name when read into ape
  
  return(mig_trees)
}