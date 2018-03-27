##### Function to compare posterior and posterior predictive summary stats #####
##### This version uses all comparison methods #####

compare_stats  <- function(run_name, descrStats, stats, randomize){ # descrStats = summary stats to compare; stats = list of dataframes of summary stats values
  
  out <- c() # create empty vector for results
  
  if ("GSI" %in% descrStats){
    out$gsi$ttest <- post.diff.calculator.ttest("GSI", stats$post$gsi, stats$post_pred$gsi, randomize) # calculate mean, sd, and t-test of difference between posterior and posterior predictive gsi values
    out$gsi$wilcoxon <- post.diff.calculator.wilcoxon("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$kstest <- post.diff.calculator.kstest("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$quantiles <- post.diff.calculator.quantiles("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$pvalue <- post.diff.calculator.pvalue("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$containment <- post.diff.calculator.containment("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$error <- post.diff.calculator.error("GSI", stats$post$gsi, stats$post_pred$gsi, randomize)
    out$gsi$total <- plyr::rbind.fill(out$gsi$ttest, out$gsi$wilcoxon, out$gsi$kstest, out$gsi$quantiles, out$gsi$pvalue, out$gsi$containment, out$gsi$error)
  }
  
  if ("SMS" %in% descrStats){
    out$sms$ttest <- post.diff.calculator.ttest("SMS", stats$post$sms, stats$post_pred$sms, randomize) # calculate mean, sd, and t-test of difference between posterior and posterior predictive sms values
    out$sms$wilcoxon <- post.diff.calculator.wilcoxon("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$kstest <- post.diff.calculator.kstest("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$quantiles <- post.diff.calculator.quantiles("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$pvalue <- post.diff.calculator.pvalue("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$containment <- post.diff.calculator.containment("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$error <- post.diff.calculator.error("SMS", stats$post$sms, stats$post_pred$sms, randomize)
    out$sms$total <- plyr::rbind.fill(out$sms$ttest, out$sms$wilcoxon, out$sms$kstest, out$sms$quantiles, out$sms$pvalue, out$sms$containment, out$sms$error)
  }
  
  if ("DIV" %in% descrStats){
    out$NDW$ttest <- post.diff.calculator.ttest("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$wilcoxon <- post.diff.calculator.wilcoxon("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$kstest <- post.diff.calculator.kstest("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$quantiles <- post.diff.calculator.quantiles("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$pvalue <- post.diff.calculator.pvalue("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$containment <- post.diff.calculator.containment("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$error <- post.diff.calculator.error("NDW", stats$post$pi_within, stats$post_pred$pi_within, randomize)
    out$NDW$total <- plyr::rbind.fill(out$NDW$ttest, out$NDW$wilcoxon, out$NDW$kstest, out$NDW$quantiles, out$NDW$pvalue, out$NDW$containment, out$NDW$error)
    
    out$NDB$ttest <- post.diff.calculator.ttest("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$wilcoxon <- post.diff.calculator.wilcoxon("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$kstest <- post.diff.calculator.kstest("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$quantiles <- post.diff.calculator.quantiles("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$pvalue <- post.diff.calculator.pvalue("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$containment <- post.diff.calculator.containment("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$error <- post.diff.calculator.error("NDB", stats$post$pi_between, stats$post_pred$pi_between, randomize)
    out$NDB$total <- plyr::rbind.fill(out$NDB$ttest, out$NDB$wilcoxon, out$NDB$kstest, out$NDB$quantiles, out$NDB$pvalue, out$NDB$containment, out$NDB$error)
    
    out$NDR$ttest <- post.diff.calculator.ttest("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$wilcoxon <- post.diff.calculator.wilcoxon("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$kstest <- post.diff.calculator.kstest("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$quantiles <- post.diff.calculator.quantiles("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$pvalue <- post.diff.calculator.pvalue("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$containment <- post.diff.calculator.containment("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$error <- post.diff.calculator.error("NDR", stats$post$pi_ratio, stats$post_pred$pi_ratio, randomize)
    out$NDR$total <- plyr::rbind.fill(out$NDR$ttest, out$NDR$wilcoxon, out$NDR$kstest, out$NDR$quantiles, out$NDR$pvalue, out$NDR$containment, out$NDR$error)
  }
  
  if ("FST" %in% descrStats){
    out$FSTP$nuc$ttest <- post.diff.calculator.ttest("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$wilcoxon <- post.diff.calculator.wilcoxon("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$kstest <- post.diff.calculator.kstest("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$quantiles <- post.diff.calculator.quantiles("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$pvalue <- post.diff.calculator.pvalue("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$containment <- post.diff.calculator.containment("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$error <- post.diff.calculator.error("FSTPN", stats$post$fst_pairwise$nuc, stats$post_pred$fst_pairwise$nuc, randomize)
    out$FSTP$nuc$total <- plyr::rbind.fill(out$FSTP$nuc$ttest, out$FSTP$nuc$wilcoxon, out$FSTP$nuc$kstest, out$FSTP$nuc$quantiles, out$FSTP$nuc$pvalue, out$FSTP$nuc$containment, out$FSTP$nuc$error)
    
    out$FSTR$nuc$ttest <- post.diff.calculator.ttest("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$wilcoxon <- post.diff.calculator.wilcoxon("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$kstest <- post.diff.calculator.kstest("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$quantiles <- post.diff.calculator.quantiles("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$pvalue <- post.diff.calculator.pvalue("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$containment <- post.diff.calculator.containment("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$error <- post.diff.calculator.error("FSTRN", stats$post$fst_range$nuc, stats$post_pred$fst_range$nuc, randomize)
    out$FSTR$nuc$total <- plyr::rbind.fill(out$FSTR$nuc$ttest, out$FSTR$nuc$wilcoxon, out$FSTR$nuc$kstest, out$FSTR$nuc$quantiles, out$FSTR$nuc$pvalue, out$FSTR$nuc$containment, out$FSTR$nuc$error)
    
    out$FSTG$nuc$ttest <- post.diff.calculator.ttest("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$wilcoxon <- post.diff.calculator.wilcoxon("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$kstest <- post.diff.calculator.kstest("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$quantiles <- post.diff.calculator.quantiles("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$pvalue <- post.diff.calculator.pvalue("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$containment <- post.diff.calculator.containment("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$error <- post.diff.calculator.error("FSTGN", stats$post$fst_global$nuc, stats$post_pred$fst_global$nuc, randomize)
    out$FSTG$nuc$total <- plyr::rbind.fill(out$FSTG$nuc$ttest, out$FSTG$nuc$wilcoxon, out$FSTG$nuc$kstest, out$FSTG$nuc$quantiles, out$FSTG$nuc$pvalue, out$FSTG$nuc$containment, out$FSTG$nuc$error)
    
    out$FSTP$hap$ttest <- post.diff.calculator.ttest("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$wilcoxon <- post.diff.calculator.wilcoxon("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$kstest <- post.diff.calculator.kstest("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$quantiles <- post.diff.calculator.quantiles("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$pvalue <- post.diff.calculator.pvalue("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$containment <- post.diff.calculator.containment("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$error <- post.diff.calculator.error("FSTPH", stats$post$fst_pairwise$hap, stats$post_pred$fst_pairwise$hap, randomize)
    out$FSTP$hap$total <- plyr::rbind.fill(out$FSTP$hap$ttest, out$FSTP$hap$wilcoxon, out$FSTP$hap$kstest, out$FSTP$hap$quantiles, out$FSTP$hap$pvalue, out$FSTP$hap$containment, out$FSTP$hap$error)
    
    out$FSTR$hap$ttest <- post.diff.calculator.ttest("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$wilcoxon <- post.diff.calculator.wilcoxon("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$kstest <- post.diff.calculator.kstest("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$quantiles <- post.diff.calculator.quantiles("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$pvalue <- post.diff.calculator.pvalue("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$containment <- post.diff.calculator.containment("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$error <- post.diff.calculator.error("FSTRH", stats$post$fst_range$hap, stats$post_pred$fst_range$hap, randomize)
    out$FSTR$hap$total <- plyr::rbind.fill(out$FSTR$hap$ttest, out$FSTR$hap$wilcoxon, out$FSTR$hap$kstest, out$FSTR$hap$quantiles, out$FSTR$hap$pvalue, out$FSTR$hap$containment, out$FSTR$hap$error)
    
    out$FSTG$hap$ttest <- post.diff.calculator.ttest("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$wilcoxon <- post.diff.calculator.wilcoxon("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$kstest <- post.diff.calculator.kstest("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$quantiles <- post.diff.calculator.quantiles("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$pvalue <- post.diff.calculator.pvalue("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$containment <- post.diff.calculator.containment("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$error <- post.diff.calculator.error("FSTGH", stats$post$fst_global$hap, stats$post_pred$fst_global$hap, randomize)
    out$FSTG$hap$total <- plyr::rbind.fill(out$FSTG$hap$ttest, out$FSTG$hap$wilcoxon, out$FSTG$hap$kstest, out$FSTG$hap$quantiles, out$FSTG$hap$pvalue, out$FSTG$hap$containment, out$FSTG$hap$error)
  }
  
  out_total <- rbind(out$gsi$total, out$sms$total, out$NDW$total, out$NDB$total, out$NDR$total, out$FSTP$nuc$total, out$FSTR$nuc$total, out$FSTG$nuc$total, out$FSTP$hap$total, out$FSTR$hap$total, out$FSTG$hap$total)
  comparison_list <- list("ttest", "wilcoxon", "kstest", "quantiles", "pvalue", "containment", "error")
  comparison_list <- unlist(lapply(comparison_list, rep, times = 11))
  stats_list <- list("GSI", "SMS", "NDW", "NDB", "NDR", "FSTPN", "FSTRN", "FSTGN", "FSTPH", "FSTRH", "FSTGH")
  stats_list <- unlist(lapply(stats_list, rep, times = 77))
  out_total <- as.data.frame(cbind(stats_list,comparison_list, out_total)) # add columns for comparison method and summary stat
  colnames(out_total) <- c("Stat", "Comp", "Locus", "Mean", "SD", "Sig")
  out_total <- apply(out_total, 2, as.character)
  
  if (randomize == TRUE){
    comp <- "random"
  }
  else {
    comp <- "same"
  }
  
  csv_name <- sprintf("%s_out_%s.csv", run_name, comp)
  write.csv(out_total, csv_name)
  
  return(out)
}