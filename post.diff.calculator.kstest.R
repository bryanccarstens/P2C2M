##### Function for calculating differences between the posterior distribution and the posterior predictive distribution using the Kolmogorov Smirnov test #####

post.diff.calculator.kstest <- function(summary_stat, post_df, post_pred_df, randomize) { # summary_stat = the statistics being used to compare the distributions, post_df = df of posterior summary stat values with loci as columns, post_pred_df = df of posterior predictive summary stat values with loci as columns

### Compare predictive data to posterior data randomly or each predictive set to the posterior set it was simulated from ###
  if (randomize == TRUE){
    post_pred_df <- dplyr::sample_n(post_pred_df, size = nrow(post_pred_df))
  }
  
#difference between posterior distribution and 1 simulation from that generation
### Each locus ###
  
  diff_df <- as.data.frame(do.call(cbind, lapply(seq(1, ncol(post_df), 1), function(x) post_pred_df[,x] - post_df[,x])))
  ks_test <- suppressWarnings(mapply(stats::ks.test, post_df, post_pred_df))
  ks_test <- ks_test[2,]
  
  diff_mean <- lapply(diff_df, mean) # calculate mean of differences for each locus
  diff_sd <- lapply(diff_df, sd) # calculate standard deviation of differences for each locus
  
  summary_df <- as.data.frame(matrix(nrow = ncol(diff_df), ncol = 0)) # create empty df for summary values
  for (z in 1:ncol(diff_df)){ # for each locus
    
    #hist.writer(summary_stat, z, diff_df[,z]) # write histogram of differences to pdf
    
    locus_summary <- list(z, unlist(diff_mean[z]), unlist(diff_sd[z]), unlist(ks_test[z])) # create list of summary values for locus
    summary_df <- rbind(summary_df, locus_summary) # add locus summary values to summary df
  }
  
  colnames(summary_df) <- c("Locus", "Mean", "SD", "Sig")
  
### Across Loci ###
  diff_mean_all <- mean(unlist(diff_df)) # calculate mean of all loci
  diff_sd_all <- sd(unlist(diff_df)) # calculate sd of all loci
  ks_test_all <- suppressWarnings(stats::ks.test(unlist(post_df), unlist(post_pred_df))) # perform t test of all differences
  ks_test_pvalue_all <- ks_test_all[2] # extract p value from t test results
  across_summary <- list("Across", unlist(diff_mean_all), unlist(diff_sd_all), unlist(ks_test_pvalue_all)) # create list of summary values across loci
  
  #hist.writer(summary_stat, "acrossloci", unlist(diff_df)) # write histogram of differences to pdf
  
  summary_df_all <- rbind(summary_df, across_summary) # add across-locus summary to summary df
  #write.csv(summary_df_all, paste(summary_stat, "_kstest_1to1_same.csv", sep = ""))
  
  return(as.data.frame(summary_df_all)) # does not include values calculated over all loci yet
}