##### Function to run P2C2M for all runs in a directory #####

p2c2m.migrate.all <- function(descrStats = c("GSI", "SMS", "DIV", "FST"), singleAllele = c(), num.reps = 1, num.sims = 1000){
  
  bayes_files <- list.files(pattern = "bayesall*")
  infiles <- list.files(pattern = "migrate_model*")
  treefiles <- list.files(pattern = "migrate_treefile*")
  
  for (q in 1:length(bayes_files)){
    p2c2m.migrate(mig_bayesall = bayes_files[q], mig_in = infiles[q], mig_tree = treefiles[q], descrStats = descrStats, singleAllele = singleAllele, num.reps = num.reps, num.sims = num.sims)
  }
}