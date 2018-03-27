##### Function for calculating the percentage of post reps that fall outside the 95% post_pred density for a given summary stat #####

pp_containment <- function(post_values, post_pred_values){ # post_df = posterior summary stat values, post_pred_df = posterior predictive summary stat values
  pred_quantiles <- stats::quantile(post_pred_values, probs = c(0.025, 0.975)) # calculate quantiles for pred values
  in_value = 0
  out_value = 0
    
  for (t in post_values){ # for each gene tree
    if (t >= pred_quantiles[1] && t <= pred_quantiles[2]){ # if post value falls within 95% of pred distribution
      in_value = in_value + 1
    }
    
    else {
      out_value = out_value + 1
    }
  }
  
  p_containment <- out_value / (in_value + out_value) # calculate percentage of post values that fall outside 95% of pred distribution
  return(p_containment)
}