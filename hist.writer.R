##### Function to write histograms to a pdf #####

hist.writer <- function(summary_statistic, locus_id, diff_data){ # summary_statistic = the stat being used to compare the distributions, locus_number = which locus to create pdf for, diff_data = data to make histogram of
  
  grDevices::pdf(file=sprintf("%s_diff_1to1_same_locus_%s.pdf",summary_statistic, locus_id)) # create pdf 
  lim <- (max(abs(diff_data))) + (0.1 * max(abs(diff_data)))
  hist(diff_data, breaks = 30, main = "Difference Between Posterior and Posterior Predictive Distributions", xlab = "", xlim = c(-lim, lim)) # write histogram of differences to pdf
  grDevices::dev.off() # close connection to pdf
  
}