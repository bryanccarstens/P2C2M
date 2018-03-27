##### Function for calcuating Bayesian posterior predictive p-values #####

pp_pvalue <- function(post_values, post_pred_values){ # post_value = value of posterior summary stat, either a single value calculated from the empirical dataset or a mean of posterior values; post_pred_values = set of summary stat values to be compared to posterior, either values for a single locus or values for all loci
  post_mean <- mean(post_values)
  u_value = 0
  l_value = 0
  
  for (val in post_pred_values){ # for each stat value
    if (val > post_mean){ # if the predictive value is greater than the posterior value
      u_value = u_value + 1 
    }
    if (val < post_mean){ # if the predictive value is less than the posterior value
      l_value = l_value + 1
    }
  }
  
  equal <- length(which(post_pred_values == post_mean)) / 2 # calculate number of predictive values that equal the posterior value
  
  #pp_pvalue <- 2 * min(u_value, l_value) / length(post_pred_values) # 2-tailed p-value
  pp_pvalue <- (u_value + equal) / length(post_pred_values)
  return(pp_pvalue)
}