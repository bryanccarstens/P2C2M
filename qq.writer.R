##### Function to write qq plots to a pdf #####

qq.writer <- function(summary_statistic, locus_id, post_df, post_pred_df){ # summary_statistic = the stat being used to compare the distributions, locus_number = which locus to create pdf for, post_df = df of posterior values, post_pred_df = df of post pred values
  
  grDevices::pdf(file=sprintf("%s_diff_1to1_same_locus_%s.pdf",summary_statistic, locus_id)) # create pdf 
  qqplot(post_df, post_pred_df, main = "Difference Between Posterior and Posterior Predictive Distributions", xlab = "") # write qqplot to pdf
  grDevices::dev.off() # close connection to pdf
  
}