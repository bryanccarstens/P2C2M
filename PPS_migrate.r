library(dplyr)

######################################################
# GET VARIABLES FROM MIGRATE FILES FOR PREDICTIVE DIST
######################################################

#read in migrate file - default = "bayesallfile"
bayesallfile<-read.table("bayesallfile", header=T, fill=T)

#get number of loci
loci<-unique(bayesallfile$Locus)
loci_n<-length(loci)


##########################################################################
#RANDOMLY SAMPLE FROM POST
#THIS MIGHT BE SKIPPED IF WE ARE USING THE ENTIRE POST
#get number of simulations to run
#gens_run<-1000 #how many generations per loci to simulate - 1000 default?
#randomly sample some number of generations
#y<-subset(bayesallfile, Locus==1)
#burnin<-ceiling(nrow(y)*0.1)
#y<-tail(y, -burnin)
#y<-sample(unique(bayesallfile$Steps),gens_run)
###########################################################################

#number of simulations to run per generation in the posterior
sims_run<-100 # how many sims per generation - 100 default?


#get sample sizes from infile
#migrate infile default = "infile"
#produces a dataframe: pops (rows) X loci (columns)
infile<-read.table("infile", fill=T, na.strings=c("","NA"))

samples<-infile[!grepl('[A-z]',infile[,2]),]
samples<-samples[-c(1:2),]
samples <- lapply(samples, function(x) as.numeric(as.character(x)))
samples<-as.data.frame(samples)


#get theta from this dataset
theta_pops<-grep("Theta_", names(bayesallfile), value=T)
theta_pops<-sort(theta_pops)

#number of populations
pops<-length(theta_pops)

#get migration events from this dataset
migration<-grep("M_", names(bayesallfile), value=T)
migration<-sort(migration)


#function to extract each migration event
make_migrations<-function (a) {
  b<-substring(a,3,5)
  b<-sub("_", " ", b)
  b<-paste("-m", b, sep=" ")
}

migrations<-lapply(migration, make_migrations)


#simulate data for each locus for each generation
for (l in loci){
  
  #conduct simualations by locus
  df<-subset(bayesallfile, Locus==l)
  #subset from post - THIS MIGHT BE SKIPPED
  #df<-df[df$Steps %in% y,]
  
  
  for (g in 1:nrow(df)) {
    
    #get max theta for each generation
    t<-select(df[g,], theta_pops)
    t<-max(as.numeric(t[1,]))
    
    #get all migrations values for each generation
    migvalues<-as.list(select(df[g,], migration))
    mig_val<-mapply(c, migrations, migvalues, SIMPLIFY=T)
    mig_val<-paste(mig_val, collapse=" ")
    
    #get population information
    total_n<-sum(samples[,l])
    popvalues<-as.list(samples[,l])
    pops_n<-paste(popvalues, collapse=" ")
    
    #paste together ms command
    ms_input1<-paste("./ms", total_n, sims_run, "-t", t, "-I", pops, pops_n, mig_val, sep=" ")
    ms_input2<-paste("-T | tail -n +4 | grep -v // >> locus_", l, "_pps.trees",  sep="")
    ms_input<-paste(ms_input1, ms_input2, sep=" ")
    
    #check output
    #write.table(data.frame(l,ms_input), file="test_msinput",append=T, col.names=F, row.names=F, quote=F)
    
    #call ms with command
    system(ms_input)
    }  
}


##################################################
#get population associations for stat calculations
##################################################

#function to fill in columns with pop assignment
repeat_last = function(x, forward = TRUE, maxgap = Inf, na.rm = FALSE) {
  if (!forward) x = rev(x)           # reverse x twice if carrying backward
  ind = which(!is.na(x))             # get positions of nonmissing values
  if (is.na(x[1]) && !na.rm)         # if it begins with NA
    ind = c(1,ind)                 # add first pos
  rep_times = diff(                  # diffing the indices + length yields how often
    c(ind, length(x) + 1) )          # they need to be repeated
  if (maxgap < Inf) {
    exceed = rep_times - 1 > maxgap  # exceeding maxgap
    if (any(exceed)) {               # any exceed?
      ind = sort(c(ind[exceed] + 1, ind))      # add NA in gaps
      rep_times = diff(c(ind, length(x) + 1) ) # diff again
    }
  }
  x = rep(x[ind], times = rep_times) # repeat the values at these indices
  if (!forward) x = rev(x)           # second reversion
  x
}

#column with population name in infile
col<-loci_n+1

#fill in pop assignment and make dataframe with assignments
infile[,col]<-repeat_last(infile[,col])
nam<-tail(infile, -3)
nam<-nam[!complete.cases(nam),]
nam<-nam[,c(1,col)]
assign<-unique(nam)

#get lists of pop assignments
pop_assign<-split(assign, assign[2])








