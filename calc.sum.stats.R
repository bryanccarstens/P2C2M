##### Function for calculating summary statistics for posterior and posterior predictive samples #####

calc.sum.stats <- function(ppt_infsites, empirical_popgenome, pred_trees, post_trees, post_samples, spal, loci, descrStats){

### Create vector and dataframes for summary stats output ###
  
  stats <- c()
  stats$post_pred$gsi <- as.data.frame(matrix(nrow = length(post_samples), ncol = 0)) # create post pred dataframe for gsi
  stats$post$gsi <- as.data.frame(matrix(nrow = length(post_samples), ncol = 0)) # create post dataframe for gsi
  stats$post_pred$sms <- as.data.frame(matrix(nrow = length(post_samples), ncol = 0)) # create post pred df for sms
  stats$post$sms <- as.data.frame(matrix(nrow = length(post_samples), ncol = 0)) # created post df for sms
 
### Calculate tree statistics ###
  
  cat("\n")
  for (locus in 1:length(loci)) { # for each locus
    print(sprintf("Analyzing Locus%s Gene Trees", locus))
    
    if ("GSI" %in% descrStats){
      stats$post_pred$gsi <- cbind(stats$post_pred$gsi, gsi.helper(spal$pred, spal$pred_popset_nosingles, pred_trees[[locus]])) # calculate gsi for predictive trees and add to stats vector
      stats$post$gsi <- cbind(stats$post$gsi, gsi.helper(spal$post, spal$post_popset_nosingles, post_trees[[locus]])) # calculate gsi for post trees and add to stats vector
    }
    
    if ("SMS" %in% descrStats){
      stats$post_pred$sms <- cbind(stats$post_pred$sms, sms.helper(spal$pred, pred_trees[[locus]])) # calculate sms for predictive trees and add to stats vector
      stats$post$sms <- cbind(stats$post$sms, sms.helper(spal$post, post_trees[[locus]])) # calculate sms for post trees and add to stats vector
    }
  }  
  
  stats$post_pred$gsi <- as.data.frame(stats$post_pred$gsi)
  stats$post$gsi <- as.data.frame(stats$post$gsi)
  stats$post_pred$sms <- as.data.frame(stats$post_pred$sms)
  stats$post$sms <- as.data.frame(stats$post$sms)
  
### Calculate diversity and FST statistics ###
  
  if ("DIV" %in% descrStats | "FST" %in% descrStats){
    cat("\nCalculating diversity/Fst statistics")
    ppt_infsites <- lapply(seq(1, length(loci), 1), function(x) PopGenome::F_ST.stats(ppt_infsites[[x]]))
    empirical_popgenome <- PopGenome::F_ST.stats(empirical_popgenome)
    
    ### Note: 0.0000000001 is added to each value in non-ratio datasets to avoid dividing by 0 ###
    if("DIV" %in% descrStats){
      stats$post_pred$pi_within <- as.data.frame((do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@nuc.diversity.within), apply, 1, mean, na.rm = TRUE)) + 0.0000000001))
      stats$post$pi_within <- as.data.frame((do.call(cbind, lapply(data.table::transpose(as.data.frame(empirical_popgenome@nuc.diversity.within)), mean, na.rm = TRUE)) + 0.0000000001))
      stats$post_pred$pi_between <- as.data.frame((do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@nuc.diversity.between), apply, 2, mean, na.rm = TRUE)) + 0.0000000001))
      stats$post$pi_between <- as.data.frame(do.call(cbind, lapply(as.data.frame(empirical_popgenome@nuc.diversity.between), mean, na.rm = TRUE)) + 0.0000000001)
      stats$post_pred$pi_ratio <- as.data.frame(stats$post_pred$pi_within / stats$post_pred$pi_between)
      stats$post$pi_ratio <- as.data.frame(stats$post$pi_within / stats$post$pi_between)
    }
    
    if("FST" %in% descrStats){
      stats$post_pred$fst_pairwise$nuc <- as.data.frame(do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@nuc.F_ST.pairwise), apply, 2, mean, na.rm = TRUE)))
      stats$post$fst_pairwise$nuc <- as.data.frame(do.call(cbind, lapply(as.data.frame(empirical_popgenome@nuc.F_ST.pairwise), mean, na.rm = TRUE)))
      stats$post_pred$fst_range$nuc <- as.data.frame(do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@nuc.F_ST.pairwise), apply, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
      stats$post$fst_range$nuc <- as.data.frame(do.call(cbind, lapply(as.data.frame(empirical_popgenome@nuc.F_ST.pairwise), function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
      stats$post_pred$fst_global$nuc <- as.data.frame(do.call(cbind, lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@nucleotide.F_ST)))
      stats$post$fst_global$nuc <- as.data.frame(transpose(as.data.frame(empirical_popgenome@nucleotide.F_ST)))
      stats$post_pred$fst_pairwise$hap <- as.data.frame(do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@hap.F_ST.pairwise), apply, 2, mean, na.rm = TRUE)))
      stats$post$fst_pairwise$hap <- as.data.frame(do.call(cbind, lapply(as.data.frame(empirical_popgenome@hap.F_ST.pairwise), mean, na.rm = TRUE)))
      stats$post_pred$fst_range$hap <- as.data.frame(do.call(cbind, lapply(lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@hap.F_ST.pairwise), apply, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
      stats$post$fst_range$hap <- as.data.frame(do.call(cbind, lapply(as.data.frame(empirical_popgenome@hap.F_ST.pairwise), function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
      stats$post_pred$fst_global$hap <- as.data.frame(do.call(cbind, lapply(seq(1, length(loci), 1), function(x) ppt_infsites[[x]]@haplotype.F_ST)))
      stats$post$fst_global$hap <- as.data.frame(data.table::transpose(as.data.frame(empirical_popgenome@haplotype.F_ST)))
    }
  }
  
  return(stats)
}


