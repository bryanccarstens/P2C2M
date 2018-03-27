##### Function to simulate posterior predictive datasets #####

pred_simulator <- function(sub_bayesallfile, pops, samples, loci, segsites_per_locus, migration, migrations, sims_run){

  for (l in loci){
    
    print(paste0("Simulating Locus", l))
    #conduct simulations by locus
    df<-subset(sub_bayesallfile, Locus==l)
    
    # get number segregating sites for locus
    locus_segsites <- segsites_per_locus[l]
    
    # start pred tree file
    treefile_name <- sprintf("Locus%s_pps.trees", l)
    pptrees_file <- file(treefile_name, "wt")
    ms_line <- paste("ms", sum(samples[,l]), nrow(df) * sims_run, "\n", sep = " ")
    writeLines(ms_line, pptrees_file)
    close(pptrees_file)
    
    # ms simulations
    for (g in 1:nrow(df)) {
      
      #get all migrations values for each generation
      migvalues<-as.list(dplyr::select_(df[g,], .dots = migration))
      mig_val<-mapply(c, migrations, migvalues, SIMPLIFY=T)
      mig_val<-paste(mig_val, collapse=" ")
      
      #get population information
      total_n<-sum(samples[,l])
      popvalues<-as.list(samples[,l])
      pops_n<-paste(popvalues, collapse=" ")
      
      #paste together ms command
      ms_input1<-paste("./ms", total_n, sims_run, "-s", locus_segsites, "-I", pops, pops_n, mig_val, sep=" ")
      ms_input2<-paste("-T | tail +4 >> Locus", l, "_pps.trees",  sep="")
      ms_input<-paste(ms_input1, ms_input2, sep=" ")
      
      #check output
      #write.table(data.frame(l,ms_input), file="test_msinput",append=T, col.names=F, row.names=F, quote=F)
      
      #call ms with command
      system(ms_input)
    }  
  }
}
