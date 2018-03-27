##### Function to get metadata and conduct posterior predictive simulations #####

pps_migrate <- function(mig_bayesall, mig_in, num.reps, num.sims, run_name){ # num.reps = the number of simulation replicates to perform for each posterior sample

######################################################
# GET VARIABLES FROM MIGRATE FILES FOR PREDICTIVE DIST
######################################################

#read in migrate file - default = "bayesallfile"
  bayesallfile<-utils::read.table(mig_bayesall, header=T, fill=T)

#get number of loci
  loci<-unique(bayesallfile$Locus)
  loci_n<-length(loci)

#number of simulations to run per generation in the posterior
  sims_run<-num.reps # how many sims per generation - 100 default?

##########################################################################
#RANDOMLY SAMPLE FROM POST
#THIS MIGHT BE SKIPPED IF WE ARE USING THE ENTIRE POST
#get number of simulations to run
  gens_sim<-num.sims #how many generations per loci to simulate - 1000 default?
#randomly sample some number of generations
  gens_run <- nrow(bayesallfile[bayesallfile$Locus == 1,]) # get number of samples in posterior
  burnin <- ceiling(gens_run * 0.1)
  non_burnin <- seq(burnin, gens_run, 1)
  if (gens_sim > length(non_burnin)){ # if number to sample > total number of samples
    post_samples <- non_burnin # simulate total number of samples
  } else{
    post_samples <- sort(sample(non_burnin, gens_sim))
  }
  post_samples_sim <- post_samples
  for (v in 1:(loci_n - 1)){
    post_samples_sim <- c(post_samples_sim, (post_samples_sim + gens_run))
  }
  post_samples_sim <- unique(post_samples_sim)
  sub_bayesallfile <- bayesallfile[post_samples_sim,]
###########################################################################

#get sample sizes from infile
#migrate infile default = "infile"
#produces a dataframe: pops (rows) X loci (columns)
  infile<-utils::read.table(mig_in, fill=T, na.strings=c("","NA"))

  samples<-infile[!grepl('[A-z]',infile[,2]),]
  samples<-samples[-c(1:2),]
  samples <- lapply(samples, function(x) as.numeric(as.character(x)))
  samples<-as.data.frame(samples)
  
# create association matrix
  pop_assign <- create_association_matrix(infile, loci_n)

#get theta from this dataset
  theta_pops<-grep("Theta_", names(sub_bayesallfile), value=T)
  theta_pops<-sort(theta_pops)

#number of populations
  pops<-length(theta_pops)

#get migration events from this dataset
  migration<-grep("xNm_", names(sub_bayesallfile), value=T)
  migration<-sort(migration)


#function to extract each migration event
  make_migrations<-function (a) {
    b<-substring(a,5,7)
    b<-sub("_", " ", b)
    b<-paste("-m", b, sep=" ")
  }
  
  migrations<-lapply(migration, make_migrations)

# Get # segregating sites for each locus from fasta files
  wd <- paste0(getwd(), "/fasta")
  wd <- paste0(wd, "_", run_name)
  empirical_popgenome <- PopGenome::readData(wd)
  segsites_per_locus <- lapply(empirical_popgenome@region.data@biallelic.sites, length)
  
# Define populations in empirical popGenome element
  unique_pops <- as.character(unique(do.call("rbind", pop_assign)[,2])) # get unique populations
  inds_list <- lapply(seq(1, length(pop_assign), 1), function(x) as.character(pop_assign[[x]][,1])) # create list of vectors containing inds for each pop
  empirical_popgenome <- PopGenome::set.populations(empirical_popgenome, inds_list) 
  cat("\n\n")
  
#simulate data for each locus for each generation
  pred_simulator(sub_bayesallfile, pops, samples, loci, segsites_per_locus, migration, migrations, sims_run)

  return(list(bayesallfile, loci, samples, pops, pop_assign, post_samples, empirical_popgenome))
}


