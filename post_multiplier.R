##### Function to repeat rows in the posterior summary statistics to match the posterior predictive simulations #####

post_multiplier <- function (stats, sims_run){ # stats = summary stats list; sims_run = number of simulations run for each posterior sample
  
  stats$post$gsi <- stats$post$gsi[rep(seq_len(nrow(stats$post$gsi)), each = sims_run),]
  stats$post$sms <- stats$post$sms[rep(seq_len(nrow(stats$post$sms)), each = sims_run),]
  
  stats$post$pi_within <- stats$post$pi_within[rep(seq_len(nrow(stats$post$pi_within)), each = sims_run),]
  stats$post$pi_between <- stats$post$pi_between[rep(seq_len(nrow(stats$post$pi_between)), each = sims_run),]
  stats$post$pi_ratio <- stats$post$pi_ratio[rep(seq_len(nrow(stats$post$pi_ratio)), each = sims_run),]
  
  stats$post$fst_pairwise$nuc <- stats$post$fst_pairwise$nuc[rep(seq_len(nrow(stats$post$fst_pairwise$nuc)), each = sims_run),]
  stats$post$fst_range$nuc <- stats$post$fst_range$nuc[rep(seq_len(nrow(stats$post$fst_range$nuc)), each = sims_run),]
  stats$post$fst_global$nuc <- stats$post$fst_global$nuc[rep(seq_len(nrow(stats$post$fst_global$nuc)), each = sims_run),]
  stats$post$fst_pairwise$hap <- stats$post$fst_pairwise$hap[rep(seq_len(nrow(stats$post$fst_pairwise$hap)), each = sims_run),]
  stats$post$fst_range$hap <- stats$post$fst_range$hap[rep(seq_len(nrow(stats$post$fst_range$hap)), each = sims_run),]
  stats$post$fst_global$hap <- stats$post$fst_global$hap[rep(seq_len(nrow(stats$post$fst_global$hap)), each = sims_run),]
  
  return(stats)
  
}