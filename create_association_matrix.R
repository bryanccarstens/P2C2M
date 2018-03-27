##################################################
#get population associations for stat calculations
##################################################

create_association_matrix <- function(infile, loci_n){
  
  #column with population name in infile
  col<-loci_n+1
  
  #fill in pop assignment and make dataframe with assignments
  infile[,col]<-repeat_last(infile[,col])
  nam<-utils::tail(infile, -3)
  nam<-nam[!stats::complete.cases(nam),]
  nam<-nam[,c(1,col)]
  assign<-unique(nam)
  
  #get lists of pop assignments
  pop_assign<-split(assign, assign[2])

  return(pop_assign)
}