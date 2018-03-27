##### Function for calculating differences between the posterior distribution and the posterior predictive distribution by determining if the differences overlap 0 #####

post.diff.calculator.quantiles <- function(summary_stat, post_df, post_pred_df, randomize) { # summary_stat = the statistics being used to compare the distributions, post_df = df of posterior summary stat values with loci as columns, post_pred_df = df of posterior predictive summary stat values with loci as columns
  
### Compare predictive data to posterior data randomly or each predictive set to the posterior set it was simulated from ###
  if (randomize == TRUE){
    post_pred_df <- dplyr::sample_n(post_pred_df, size = nrow(post_pred_df))
  }
  
#difference between posterior distribution and 1 simulation from that generation
### Each locus ###
  diff_df <- as.data.frame(do.call(cbind, lapply(seq(1, ncol(post_df), 1), function(x) post_pred_df[,x] - post_df[,x])))
  
  diff_mean <- lapply(diff_df, mean) # calculate mean of differences for each locus
  diff_sd <- lapply(diff_df, sd) # calculate standard deviation of differences for each locus
  
  quant_list <- lapply(diff_df, stats::quantile, probs = c(0.025, 0.975))
  quant_sig <- list()
  
  for (q in 1:length(quant_list)){
    if (quant_list[q][[1]][1] < 0 & quant_list[q][[1]][2] > 0){
      quant_sig[q] <- 1
    }
    
    else {
      quant_sig[q] <- 0
    }
  }
  
  
  summary_df <- cbind(seq(1, ncol(post_df)), diff_mean, diff_sd, quant_sig)
  
  
  colnames(summary_df) <- c("Locus", "Mean", "SD", "Sig")
  
### Across Loci ###
  diff_mean_all <- mean(unlist(diff_df)) # calculate mean of all loci
  diff_sd_all <- sd(unlist(diff_df)) # calculate sd of all loci
  
  across_quant <- stats::quantile(unlist(diff_df), probs = c(0.025, 0.975))
  
  if (across_quant[1] < 0 & across_quant[2] > 0){
    across_sig <- 1
  }
  
  else {
    across_sig <- 0
  }
  
  across_summary <- list("Across", unlist(diff_mean_all), unlist(diff_sd_all), across_sig) # create list of summary values across loci
  
  #hist.writer(summary_stat, "acrossloci", unlist(diff_df)) # write histogram of differences to pdf
  
  summary_df_all <- rbind(summary_df, across_summary) # add across-locus summary to summary df
  #write.csv(summary_df_all, paste(summary_stat, "_quants_1to1_same.csv", sep = ""))
  
  return(as.data.frame(summary_df_all)) # does not include values calculated over all loci yet
}