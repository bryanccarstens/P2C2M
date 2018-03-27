##### Function to run P2C2M for Migrate-n #####

p2c2m.migrate <- function(mig_bayesall = "bayesallfile.gz", mig_in = "infile", mig_tree = "treefile", descrStats = c("GSI", "SMS", "DIV", "FST"), singleAllele = c(), num.reps = 1, num.sims = 1000, print_out = FALSE){ # descrStats = the summary stats to use for comparisons; singleAllele = populations with a single individual in the dataset; num.reps = the number of simulation replicates to perform for each posterior sample; randomize = whether to randomly compare predictive datasets to posterior datasets or compare predictive datasets to the posteriors they were simulated from
  
  run_name <- paste(strsplit(mig_tree, "_")[[1]][3:5], collapse = "_") # get run name - take out in final version
  
  pps <- pps_migrate(mig_bayesall, mig_in, num.reps, num.sims, run_name) # get metadata and conduct posterior predictive simulations
  
  bayesallfile <- pps[1][[1]]
  loci <- pps[2][[1]]
  samples <- pps[3][[1]]
  pops <- pps[4][[1]]
  pop_assign <- pps[5][[1]]
  post_samples <- pps[6][[1]]
  empirical_popgenome <- pps[7][[1]]
  
  stats <- sum_stats(bayesallfile, mig_tree, loci, samples, pops, post_samples, descrStats, pop_assign, singleAllele, empirical_popgenome) # calculate summary statistics
  
  cat("\nComparing posterior and posterior predictive data")
  
  stats <- post_multiplier(stats, num.reps) # acount for simulating multiple reps from the same posterior sample
  
  out <- compare_stats(run_name, descrStats, stats, randomize = TRUE) # compare posterior and posterior predictive summary statistics
  out <- compare_stats(run_name, descrStats, stats, randomize = FALSE)
  
  cat("\n")
  
  if (print_out == TRUE){
    print(out)
  }
}