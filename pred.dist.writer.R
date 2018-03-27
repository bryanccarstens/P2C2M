##### Function to plot posterior predictive distributions #####

pred.dist.writer <- function(summary_statistic, locus_id, post_values, post_pred_values){
  
  post_mean <- mean(post_values)
  
  grDevices::pdf(file=sprintf("%s_diff_1to1_same_locus_%s.pdf",summary_statistic, locus_id)) # create pdf 
  lim <- (max(abs(post_pred_values))) + (0.1 * max(abs(post_pred_values)))
  hist(post_pred_values, breaks = 30, main = "Posterior Predictive Distribution", xlab = "", xlim = c(-lim, lim)) # write histogram of differences to pdf
  abline(v = post_mean, col = red) # add line for posterior value
  grDevices::dev.off() # close connection to pdf
  
}