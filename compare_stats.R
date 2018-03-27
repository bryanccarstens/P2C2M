##### Function to compare posterior and posterior predictive summary stats #####

compare_stats  <- function(descrStats, stats, sims_run, randomize){ # descrStats = summary stats to compare; stats = list of dataframes of summary stats values; sims_run = number of simulations run for each posterior sample; randomize = whether to randomly compare predictive datasets to posterior datasets or compare predictive datasets to the posteriors they were simulated from
  
  stats <- post_multiplier(stats, sims_run) # acount for simulating multiple reps from the same posterior sample
  
  out <- c() # create empty vector for results
  
  if ("GSI" %in% descrStats){
    out$gsi <- post.diff.calculator("GSI", stats$post$gsi, stats$post_pred$gsi, randomize) # calculate mean, sd, and t-test of difference between posterior and posterior predictive gsi values
  }
  
  if ("SMS" %in% descrStats){
    out$sms <- post.diff.calculator("SMS", stats$post$sms, stats$post_pred$sms, randomize) # calculate mean, sd, and t-test of difference between posterior and posterior predictive sms values
  }
  
  if ("DIV" %in% descrStats){
    out$NDW <- post.diff.calculator("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDB <- post.diff.calculator("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDR <- post.diff.calculator("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
  }
  
  if ("FST" %in% descrStats){
    out$FSTP$nuc <- post.diff.calculator("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTR$nuc <- post.diff.calculator("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTG$nuc <- post.diff.calculator("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTP$hap <- post.diff.calculator("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTR$hap <- post.diff.calculator("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTG$hap <- post.diff.calculator("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
  }
  
  return(out)
}