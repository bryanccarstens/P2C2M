#############################
#   CALCULATE SUM STATS   ###
#############################

#Script to get rid of extra info in migrate treefile
#read posterior dist from migrate file
migrate_trees<-read.tree("treefile")


for (loc in loci) {
  
  file=paste("locus_",loc,"_pps.trees",sep="")
  pred_trees<-read.tree(file)
  
  mig_trees<-#separate trees by locus
  
  for (g in pred_trees){
    
    #GSI script
    #put all GSI values in dataframe called "post_pred_dist")
  }
  
    for (h in mig_trees){
  
    #call GSI script
    #put all GSI values in dataframe called "post_dist")
  }
    
  ###############
  # COMPARISONS #
  ###############

  #get one sample from each generation of post_pred_dist for one to one comparison
  post_pred_dist_1 <- post_pred_dist[seq(1, nrow(post_pred_dist),sims_run),]
  post_pred_dist_1 <- as.data.frame(post_pred_dist_1)
  #difference between posterior distribution and 1 simulation from that generation
  diff = post_dist - post_pred_dist_1
  diff_mean <- mean(diff[,1])
  diff_sd <- sd(diff[,1])
  diff_test<-t.test(diff)
  
  #write histograms to files
  pdf(file=paste("diff_1to1_same_locus_",loc,".pdf",sep=""))
  hist(diff)
  dev.off()
  
  write.csv(as.data.frame(loc, diff_mean, diff_sd, diff_test$p.value), "diff_1to1_same.csv", append=T)


  #difference in posterior distribution to random simulation from post pred dist
  post_number<-nrow(post_dist)
  post_pred_dist_n <- sample_n(post_pred_dist, post_number)
  diff =  post_dist - post_pred_dist_n
  write.csv(as.data.frame(loc, diff_mean, diff_sd, diff_test$p.value), "diff_1to1_random.csv",append=T)
  
  #write histograms to files
  pdf(file=paste("diff_1to1_random_locus_",loc,".pdf",sep=""))
  hist(diff)
  dev.off()

}


