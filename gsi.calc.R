##### Function for calculating the genealogical sortng index (gsi) #####

gsi.calc <- function(tree, grp){ # tree = gene or species tree to calculate gsi within, grp = the group to calculate gsi for
  
  n <- length(grp) - 1 # calculate n, the number of elements in the group -1
  grp.index <- match(grp, tree$tip.label) # get the ape package tree indices for the selected group
  all_paths <- list() # create empty list to store node paths in

### Obtain the node path/set for each element in group ###    
  for (i in 1:length(grp.index)){ # for each element in the selected group
    
    tip_num <- grp.index[i] # record the group index
    tip_path <- list() # create empty list for the nodepath
    
    # Obtain the node path #     
    while (all(!is.na(tip_path))){ # While there are no NAs returned (NA is returned once there is no parent node)
      tip_num <- tree$edge[,1][match(tip_num, tree$edge[,2])] # get tree index of parent node
      tip_path <- append(tip_path, tip_num, after = length(tip_path)) # append parent node index to tip_path list
    }
    
    tip_path <- sort(as.numeric(tip_path[!is.na(tip_path)])) # remove NAs from the list and sort it
    
    listname <- sprintf("node_path%s", i) # create variable name for tip_path list that is particular to the current element of the group
    all_paths[[i]] <- assign(listname, tip_path) # assign variable name to tip_path list
    
  }
  
### Determine which nodes connect all the elements in the group ###
  max_int <- max(Reduce(intersect, all_paths)) # get the highest node index that is shared among all elements of group
  overlap <- unique(unlist(all_paths)) # combine unique node indices from all element node paths
  overnodes <- as.list(overlap[match(max_int, overlap):length(overlap)]) # subset list starting at list index of the highest shared node in the combined list

### Calculate the gsi ###    
  degree <- table(match(as.list(tree$edge), overnodes)) # get the degree of each node
  obs_gs <- n / (sum(degree - 2) + sum(degree == 2)) # calculate the observed gs. Note: Ape gives the root node a degree of two, but it should be counted as three degrees. This is accounted for in the "sum(degree == 2)" statement
  degree_total <- table(tree$edge)[seq(Ntip(tree) + 1, Ntip(tree) + Nnode(tree))] # get the degrees of all nodes
  min_gs <- n / (sum(degree_total - 2) + sum(degree_total == 2)) # calculate the minimum gs for the tree. Note: Ape gives the root node a degree of two, but it should be counted as three degrees. This is accounted for in the "sum(degree == 2)" statement
  gsi <- (obs_gs - min_gs) / (1 - min_gs) # calculate the gsi
  return(gsi)

}



