##### Function for calculating the standard error between the posterior and posterior predictive values #####

pp_error <- function(post_values, post_pred_values){
  post_mean <- mean(post_values)
  pred_mean <- mean(post_pred_values)
  pred_stdev <- sd(post_pred_values)
  
  pp_error <- abs(post_mean - pred_mean) / pred_stdev
  return(pp_error)
}