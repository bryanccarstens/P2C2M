##### Function for calculating Slatkin and Maddison's s #####

sms.calc <- function(tree, taxon_association){ # tree = the tree to evaluate, as an ape phylo object, taxon_association = association table for samples and populations/species
  
### get necessary tree information ###
  tip_labels <- tree$tip.label # get tip labels
  tip_index <- match(tip_labels, tree$tip.label) # get tip label indices
  node_list <- tree$edge[,1][match(tip_index, tree$edge[,2])] #  get nodes attached to tips
  tip_index_table <- cbind(tip_index, tip_labels) # create lookup table to convert ape indices to tip names
  node_index <- seq(nrow(tip_index_table) + 1, length(tree$tip.label) + tree$Nnode)
  node_table <- cbind(node_index, node_index)
  tip_index_table <- rbind(tip_index_table, node_table)
  
### set up table of node states ###
  states_table <- cbind(as.character(taxon_association[,1]), NA) # create table for storing state values
  states_table <- rbind(states_table, cbind(node_index, NA))
  grps <- taxon_association[,2][match(tree$tip.label, taxon_association[,1])] # get states of tip labels
  states_table[,2][match(tree$tip.label, taxon_association[,1])] <- as.character(grps) # change states of tip labels
  s = 0 # set state change counter to 0
  
### calculate s over each node ###
  while (length(node_list) > 1){ # while at least one node occurs
    node_counts <- plyr::count(as.factor(node_list)) # count the number of times each node appears in list (these nodes connect the previous level)
    current_nodes <- as.character(node_counts[,1][which(node_counts[,2] > 1)]) # get which nodes occur more than once
    node_list <- node_list[which(node_list %in% current_nodes == FALSE)] # remove current nodes from node list
    
    for (n in 1:length(current_nodes)){ # for each node
      node_tips <- tree$edge[,2][which(current_nodes[n] == tree$edge[,1])] # get the tips that are connected by the node
      states_1 <- states_table[,2][which(match(states_table[,1], tip_index_table[,2][match(node_tips, tip_index_table[,1])]) != "NA")][1] # get state of first tip
      states_1 <- as.character(strsplit(as.character(states_1), ",")[[1]]) # format state to separate state values
      states_2 <- states_table[,2][which(match(states_table[,1], tip_index_table[,2][match(node_tips, tip_index_table[,1])]) != "NA")][2] # get state of second tip
      states_2 <- as.character(strsplit(as.character(states_2), ",")[[1]]) # format state to separate state values
      
      if (length(intersect(states_1, states_2)) > 0){ # if a state appears more than once in the tips states
        new_state <- intersect(states_1, states_2) # get the state values that occur most often in both states
      }
      
      else { # if no state appears more than once
        new_state <- c(states_1, states_2) # combine the two states to create a new state
        s = s + 1 # add one to the state change counter
      }
      
      new_state <- gsub("\\s+", "", paste(new_state, collapse = ",")) # format new state so it can be added to state table
      states_table[,2][which(current_nodes[n] == states_table[,1])] <- new_state # change state of node to combined tip states
      
    }
    
    new_nodes <- na.omit(tree$edge[,1][match(current_nodes, tree$edge[,2])]) # get nodes connecting nodes just analyzed
    node_list <- c(node_list, new_nodes) # add new nodes to node list
  }
  
  return(s) 
  
}
