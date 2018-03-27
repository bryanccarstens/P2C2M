##### Function write_parm writes a parameter file for migrate #####
# run_name is the name of the run you would like to perform with migrate (string)
# currently, to adjust the prior distribution types (uniform, exponential, etc), you will have to adjust them in lines 522-523 rather than in the variable definitions below

def write_parm(run_name):
	### Define variables
	# 0 = run_name
	# 1
	theta_pmin = 0.000000 # min for the uniform theta prior. default = 0.000000
	# 2
	theta_pmax = 0.100000 # max for the uniform theta prior. default = 0.100000
	# 3
	theta_palpha = 0.010000 # window size for the uniform theta prior. default = 0.010000
	# 4
	m_pmin = 0.000000 # min for the migration uniform prior. default = 0.000000
	# 5
	m_pmax = 1000.000000 # max for the uniform migration prior. default = 1000.000000
	# 6
	m_palpha = 100.000000 # window size for the uniform migration prior. default = 100.000000
	# 7
	chain_inc = 100 # sampling increment of the markov chain. default = 100
	# 8
	chain_sam = 5000 # number of samples/genealogies recorded. default = 5000
	# 9
	burnin = 10000 # number of samples to discard at beginning of chain. default = 10000
	
	parm_name = 'parmfile_{0}'.format(run_name) # create name for parmfile
	parm = open(parm_name, 'w+') #create parmfile
	
	# The next commands write the entire parmfile. Note: these are broken up into segments because the .format function gets confused when it sees othe {} symbols in the string. The same happens when using %
	parm.writelines('\
################################################################################\n\
# Parmfile for Migrate 3.6.11-June-18-15 [do not remove these first TWO lines]\n\
# generated automatically on\n\
# Wed Nov  1 15:57:51 2017\n\
#\n\
# please report problems to Peter Beerli\n\
#  email: beerli@fsu.edu\n\
#  http://popgen.sc.fsu.edu/migrate.html\n\
################################################################################\n\
#\n\
################################################################################\n\
# General options\n\
################################################################################\n\
#\n\
# Interactive or batch job usage\n\
#   Syntax: menu= < YES | NO > \n\
# For batch runs it needs to be set to NO\n\
menu=YES\n\
#\n\
# Specification of length of names of indiviudals\n\
#    Syntax: nmlength=<INTEGER between 0 .. 30>\n\
nmlength=10\n\
#\n\
#\n\
################################################################################\n\
# Data options\n\
################################################################################\n\
#\n\
# Several different main datatypes are possible:\n\
# INFINITE ALLELE: usable for electrophoretic markers,\n\
#                  other markers with unknown mutation model\n\
# STEPWISE MUTATION: usable for microsatellite data or\n\
#                  other markers with stepwise change\n\
#                  from one allele to another\n\
#                  [singlestep versus multistep model, see micro-submodel option]\n\
# FINITE SITES MUTATION: standard DNA/RNA sequence mutation\n\
#                  model, usable for DNA or RNA contiguous\n\
#                  sequences or varialbe sites only (SNP)\n\
# GENEALOGY SUMMARY: reanalyzing an old migrate run\n\
#\n\
#-------------------------------------------------------------------------------\n\
# INFINITE ALLELE\n\
#  Syntax: datatype=ALLELICDATA \n\
#          include-unknown=<YES | NO> with YES unknown alleles\n\
#                are included into analysis, NO is the default\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# STEPWISE MUTATION\n\
#  Syntax: datatype=<MICROSATELLITEDATA | BROWNIANDATA\n\
#                MICRO specifies the standard stepwise mutation\n\
#                model, the BROWNIAN is an approximation to this\n\
#          micro-submodel=<1|2:{tune,pinc}>\n\
#                 1 means singlestep mutation model (this is the default and the standard\n\
#                 2 is the Multistep model (see Watkins 2007 TPB, section 4.2) it needs\n\
#                   two parameters: tune specifies how close the model is to a singlestep model\n\
#                   so tune=0 --> singlestep, tune=1 --> infinite allele model;\n\
#                   the second parameter defines the probability that the repeat number\n\
#                   is increasing, this value cannot be larger than 0.666, I suggest 0.5.\n\
#                   Example: micro-submodel=2:{0.5,0.5}\n\
#          micro-threshold=<INTEGER> Default is 10 [MICRO only, NEEDS TO BE EVEN!],\n\
#                smaller values speed up analysis, but might also\n\
#                crash, large values slow down analysis considerably.\n\
#                Change this value only when you suspect that your\n\
#                data has huge gaps in repeat length.\n\
#          include-unknown=<YES | NO> with YES unknown alleles\n\
#                are included into analysis, NO is the default\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# FINITE SITES MUTATION\n\
#  Syntax: datatype=<SEQUENCEDATA | NUCLEOTIDE | UNLINKEDSNPS | ANCESTRAL\n\
#         SEQENCEDATA: typical linked stretches of DNA, for example mtDNA\n\
#         NUCLEOTIDE: linked DNA stretches, all invariable sites removed\n\
#         UNLINKEDSNPS: each variable site is a locus, DO NOT USE THIS YET\n\
#         ANCESTRAL: instead taking into account all posible states, use\n\
#                use only the most likely state probability, DON\'T USE THIS YET\n\
#\n\
#          freqs-from-data=<YES | NO: freq(A), freq(C), freq(G), freq(T)>\n\
#                calculate the prior base frequencies from the data,\n\
#                or specify the frequencies\n\
#          ttratio=<RATIO1 RATIO2 ....> Default is 2.0,\n\
#                ratio between transitions and transversions.\n\
#          seq-error=<VALUE> Default is 0.0, typical values for ABI 3700 \n\
#                sequencers after base calling are around 0.001 (1/650)\n\
#          categories=<VALUE:CATFILE> The categories are integers or letters\n\
#                specified in file called CATFILE, this assumes that all\n\
#                sites belong to known categories, this can be used to\n\
#                weight third positions etc.\n\
#          rates=<VALUE1 VALUE2 ...> the rates are specified arbitrarily or\n\
#                then are from a Gamma distribution with alpha=x, currently\n\
#                 the alpha value gets lost and is not recorded in the parmfile\n\
#          prob-rates=<RATE2 RATE1 ... > These rates can be arbitrary or \n\
#                generated with gamma-deviated rates and then are derived\n\
#                using Laguerre\'s quadrature, this should get better\n\
#                results than equal probability methods.\n\
#          autocorrelation=<NO | YES:VALUE> Default is NO\n\
#                autocorrelation makes only sense with rates,\n\
#                VALUE should be >1.0\n\
#          weights=<NO | YES:WEIGHTFILE> The weights are specified\n\
#                in file called WEIGHTFILE, this assumes that all sites\n\
#                belong to known weights, this can be used to weight\n\
#                portions of the sequence etc.\n\
#          interleaved=<YES | NO> Use either an interleaved or \n\
#                non-interleaved format. Default is NO,\n\
#                interleaved=YES is discouraged\n\
#          fast-likelihood=<YES | NO> Default is YES, use NO when you\n\
#                have many hundred individuals and get strange errors\n\
#                during a run, NO is scaling the conditional likelihood\n\
#                so that very small values are >0.00000\n\
#          inheritance-scalars={values for each locus}\n\
#                these values are multiplied with Theta, for example having\n\
#                two autosomal and a locus on X- and one on Y-chromosome we would give \n\
#                inheritance-scalars={1 1 0.75 0.25}\n\
#                [if all loci have the same scalar, just use {1}, even for many loci]]\n\
#          population-relabel={assignment for each location in the infile}\n\
#                example is population-relabel={1 2 2}\n\
#          random-subset=number<:seed>\n\
#                allows to subset the dataset randomly, if number > sample in population\n\
#                all samples are taken, if number is smaller then the pop sample is shuffled and\n\
#                and the first number samples are taken.\n\
#                the random number seed guarantees that the\n\
#                same subset is chosen in different runs\n\
#          usertree=<NO | UPGMA | AUTOMATIC | TREE:TREEFILE | DISTANCE:DISTFILE | RANDOM>\n\
#                Default is RANDOM, NO delivers a UPGMA tree using the data\n\
#                with TREE and DISTANCE the user needs to \n\
#                give a usertreefile or a pairwise distance file, with RANDOM\n\
#                a random tree will be the starting tree\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
#\n\
datatype=SequenceData\n\
ttratio=2.000000 \n\
freqs-from-data=YES\n\
seqerror-rate=0.000000\n\
categories=1 #no categories file specified\n\
rates=1: 1.000000 \n\
prob-rates=1: 1.000000 \n\
autocorrelation=NO\n\
weights=NO\n\
interleaved=NO\n\
fast-likelihood=NO\n\
inheritance-scalars={1.00000000000000000000}\n\
population-relabel={1}\n\
usertree=RANDOMTREE\n\
#\n\
')
	parm.writelines('\
################################################################################\n\
# Input options\n\
################################################################################\n\
#\n\
# input file location\n\
#   Syntax infile=FILEPATH\n\
infile=migrate_{0}\n\
#\n\
# Random number seed specification\n\
#   Syntax random-seed=<AUTO | OWN:< seedfile | value >\n\
#      AUTO           uses computer system clock to generate seed\n\
#      OWN:seedfile   uses file seedfile with random number seed\n\
#      OWN:value      uses number value for seed\n\
random-seed=AUTO #OWN:377392319\n\
#\n\
# Specify the title of the run, will be overridden by title in datafile\n\
#    Syntax: title=title text [up to 80 characters]\n\
title=\n\
#\n\
#\n\
################################################################################\n\
# Output options\n\
################################################################################\n\
#\n\
# Progress report to the window where the program was started\n\
#    Syntax: progress=<NO | YES | VERBOSE>\n\
#          NO       nothing is printed to the console\n\
#          YES      some messages about progress are reported [default]\n\
#          VERBOSE  more messages are reported to console\n\
progress=YES\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Recording messages to screen into logfile\n\
#   Syntax logfile=<NO | YES:logfilename>\n\
#       NONE     no recording of progress\n\
#       logfilename  path to logfile\n\
logfile=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Print the data as read into the program\n\
#   Syntax print-data=<NO | YES>\n\
print-data=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Print output to file [default is outfile]\n\
#   Syntax outfile=outfilename\n\
outfile=migrate_out_{0}\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Print output to a PDF file [default is outfile.pdf]\n\
#   Syntax pdf-outfile=outfilename.pdf\n\
pdf-outfile=migrate_out_{0}.pdf\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Report M (=migration rate/mutation rate) instead of 4Nm or 2 Nm or Nm\n\
#   Syntax use-M=<NO | YES> Default is YES, the name 4Nm is ambiguous\n\
#      for non-diploid data\n\
use-M=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
'.format(run_name))
	parm.writelines('\
# Plotting parameters: migration versus population size, such that Theta1 x immigration_.1\n\
# this shows the sum of all imigrations int a population\n\
#   Syntax plot=<NO | YES:<BOTH | OUTFILE>:<LOG | STD>\n\
#          {x-start, x-end, y-start, y-end}:<N | M>:interval>\n\
#      NO   do not show a plot\n\
#      YES  show plot with following specifications\n\
#           BOTH    print raw coordinates into MATHFILE and plot to OUTFILE\n\
#           OUTFILE plot only to OUTFILE\n\
#             LOG   scaling of both axes\n\
#             STD   non-log scaling\n\
#             {...} plot range of both parameters\n\
#             N     use xNm to plot immigration, x=<1,2,3,4>\n\
#                   depending on the inheritance characteristic of the data\n\
#             M     plot migration rate/mutation rate as immigration axis\n\
#             interval the plot range is broken up into interval intervals\n\
plot=NO\n\
#\n\
# Print plot data into a file \n\
#   Syntax: mathtfile=mathfile the values are printed in a mathematica readable way\n\
mathfile=mathfile\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Profile likelihood for each estimated parameter\n\
#   Syntax profile=<NONE | <ALL | TABLES | SUMMARY>:\n\
#               <PRECISE | DISCRETE | QUICK | FAST>  >\n\
#      NONE    do not calculate profile likelihoods\n\
#      ALL     print individual profile tables and summary [default]\n\
#      TABLES  show only tables and no summary\n\
#      SUMMARY show only summary\n\
#           PRECISE  evaluate profile likelihood at percentiles [Default]\n\
#           QUICK    assumes that there is no interaction of parameters\n\
#           FAST     same as QUICK except in last calculation cycle assumes interaction\n\
#           DISCRETE uses fixed mutipliers: 0.02,0.1,0.2,0.5,1,2,5,10,50\n\
profile=ALL:PRECISE\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Print tree into treefile\n\
#   Syntax print-tree=< NONE | <ALL | BEST | LASTCHAIN:Increment>:treefile >\n\
#         NONE no tree printed [Default, and only choice using parallel\n\
#         ALL  print all visited genealogies [careful this will be huge]\n\
#         BEST print only the best tree visited\n\
#         LASTCHAIN print all trees in last chain\n\
#         with increment INCREMENT\n\
print-tree=ALL\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# write intermediate minimal statistics into a file for later use\n\
#   Syntax write-summary=<NO | YES:SUMFILE >\n\
#                Default is NO, with YES the user needs to \n\
#                give a file to record the summary statistics\n\
write-summary=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Likelihood ratio test\n\
#   Syntax l-ratio=<NO | YES:values_to_test>\n\
#       Values_to_test are compared to the values generated in the run\n\
#   values_to_test={ab..bbab..ba ... a}\n\
#        the {} is a square matrix with values for the population sizes\n\
#        on the diagonal and migration rates off-diagonal\n\
#        the values a for the diagonal can be any of these:\n\
#        number  constant, the value is for example 0.002\n\
#        *       free to vary, the default is * for every parameter\n\
#        m       mean of theta, this can be a subgroup of all thetas\n\
#                for example the theta 1-3 are averaged and thetas 4,5 are estimated\n\
#        the values b for the migration rates can be any of these:\n\
#        number  constant, the value is for example 45.0 or 0.0\n\
#        *       free to vary, the default is * for every parameter\n\
#        m       mean of M_ij, this can be a subgroup of migration rates\n\
#                for example the M_1-3i are averaged and M_4,5i are estimated\n\
#        M       means of 4Nm (diploid), 2Nm (haploid), Nm (mtDNA, Y-chromosome)\n\
#        s       symmetric migration rates M\n\
#        S       symmetric migrants 4Nm\n\
#        an example for 5 populations could look like this:\n\
#        l-ratio=YES:{*s00s s*s00 0s*s0 00s*s s00s*\n\
#        this describes a circular stepping stone model with 5 symmetric rates\n\
#         and independent sizes, a very basic stepping stone with 2 parameters would\n\
#        look like this l-ratio=YES:{mm00m mmm00 0mmm0 00mmm m00mm}\n\
#        [The L-RATIO statement can be repeated]\n\
#  Default: l-ratio=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# AIC model selection [do not use yet, will come in Summer 2004]\n\
#   Syntax aic-modeltest=<NO | YES:<FAST | EXHAUSTIVE>>\n\
#       FAST        [do not use yet]\n\
#       EXHAUSTIVE  [do not use yet]\n\
aic-modeltest=NO\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
# Print a histogram of the time of migration events for each M(i->j)\n\
#    Syntax  mig-histogram=<NO | <ALL | MIGRATIONEVENTSONLY>:binsize:mighistfile >\n\
#         NO            do not record any events\n\
#         ALL           record migration and coalescence event\n\
#         MIGRATIONEVENTSONLY record only migration events\n\
#         binsize has to be in mutation units, with an average Theta=0.01 try 0.001\n\
# Print a histogram of the parameters through time (skyline plot)\n\
#    Syntax  skyline=<NO | YES>:binsize:skylinefile >\n\
#         NO            do not calculate parameter estimates through time\n\
#         YES           calculate parameters through time\n\
#         binsize has to be in mutation units, with an average Theta=0.01 try 0.001\n\
#         If the interval is too fine the output will be very noisy\n\
mig-histogram=NO\n\
skyline=NO #needs mig-histogram=ALL:...\n\
#\n\
#\n\
################################################################################\n\
# Parameter start settings\n\
################################################################################\n\
#\n\
#   Syntax: theta=<FST | OWN:<{value} | {value1, value2, ...., valuen} | NRANDOM:{mean std} | URANDOM{min,max}>\n\
#      migrationt=<FST | OWN:<{value} | {value1, value2, ...., valuen} | NRANDOM:{mean std} | URANDOM{min,max}>\n\
#        FST     starting parameter are derived from\n\
#                an FST-like calculation (Beerli&Felsenstein 1999\n\
#        OWN     starting values are supplied by user\n\
#           {value}   if only one value is supplied then all population\n\
#                     have the same starting value\n\
#           {value1, value2, ..., valuen} each population has its\n\
#                     own starting value, if the number of values is\n\
#                     insuffient, then the last value is the template\n\
#                     for the remaining populations\n\
#        NRANDOM  starting parameter is drawn randomely from a Normal distribution\n\
#           {mean std} with mean and standard deviation\n\
#        URANDOM  starting parameter is drawn randomely from a Uniform distribution\n\
#           {min max} with minimum and maximum values\n\
theta=FST\n\
migration=FST\n\
#\n\
#-------------------------------------------------------------------------------\n\
# Mutation rate modifiers\n\
#\n\
#   Syntax: mutation=<NOGAMMA | CONSTANT | ESTIMATE | GAMMA:alpha | OWN:loci: rate1 rate2 ... rate_loci>\n\
#      NOGAMMA      all loci have same mutation rate\n\
#      CONSTANT     all loci have same mutation rate\n\
#      ESTIMATE     BAYESIAN estimate: mutation rate is drawn from prior\n\
#      GAMMA:alpha  ML estimate: mutation rate has Gamma distribution with alpha\n\
#      OWN          mutation rate is different for every locus, but fixed\n\
#         :loci: rate1, ...     number of loci, rate of locus 1, locus 2 etc.\n\
#      DATA         mutation rate modifier is deducted from loci in the data\n\
#                   using Watterson\'s Theta and then scaling all rates Theta_locus/mean(Theta_loci\n\
mutation=CONSTANT\n\
#\n\
# Treatment of inviariant sequence loci\n\
# Syntax: analyze-loci=<A | F | V>\n\
#         A = analyze all loci (Default!)\n\
#         F = analyze all variable loci and ONE invariant and extrapolate\n\
#         V = analyze only variable loci\n\
#analyze-loci=A\n\
#\n\
#-------------------------------------------------------------------------------\n\
# FST model\n\
#\n\
fst-type=THETA\n\
#\n\
#-------------------------------------------------------------------------------\n\
# Custom migration model\n\
#\n\
#    Syntax: custom-migration={ab..bbab..ba ... a}\n\
#        the {} is a square matrix with values for the population sizes\n\
#        on the diagonal and migration rates off-diagonal\n\
#        the values a for the diagonal can be any of these:\n\
#        c       constant, the value needs to be defined in the theta option\n\
#        *       free to vary, the default is * for every parameter\n\
#        m       mean of theta, this can be a subgroup of all thetas\n\
#                for example the theta 1-3 are averaged and thetas 4,5 are estimated\n\
#        the values b for the migration rates can be any of these:\n\
#        c       constant, the value needs to be defined in the migration option\n\
#        *       free to vary, the default is * for every parameter\n\
#        m       mean of M_ij, this can be a subgroup of migration rates\n\
#                for example the M_1-3i are averaged and M_4,5i are estimated\n\
#        M       means of 4Nm (diploid), 2Nm (haploid), Nm (mtDNA, Y-chromosome)\n\
#        s       symmetric migration rates M\n\
#        S       symmetric migrants 4Nm\n\
#        an example for 5 populations could look like this:\n\
#        custom-migration={*s00s s*s00 0s*s0 00s*s s00s*\n\
#        this describes a circular stepping stone model with 5 symmetric rates\n\
#         and independent sizes, a very basic stepping stone with 2 parameters would\n\
#        look like this custom-migration={mm00m mmm00 0mmm0 00mmm m00mm}\n\
custom-migration={**}\n\
#\n\
# Influence of geography on migration rate\n\
# a distance matrix between populations changes the migration rate matrix so that\n\
# (genetic?) migration rates =  inferred migration rate / distance ~ a dispersion coefficient\n\
# the geofile contains a number of populations, names for populations (10 characters), they\n\
# need to be in order of the dataset. And the distances between the populations, they do not\n\
# need to be symmetric\n\
#    Syntax: geo:<NO | YES:filename>\n\
#             NO       distances among populations are considered to be 1 [all equal]\n\
#             YES      distances are read from a file\n\
geo=NO\n\
#\n\
#\n\
')
	parm.writelines('\
################################################################################\n\
# Search strategies\n\
################################################################################\n\
#\n\
# MCMC Strategy method\n\
#    Syntax: bayes-update=< NO | YES>\n\
#        NO      maximum likelihood method\n\
#        YES     Bayesian method\n\
# Some of the options are only available in one or other mode\n\
# BAYESIAN OPTIONS\n\
#        bayes-updatefreq=VALUE \n\
#            VALUE      is a ratio between 0 and 1\n\
#                       ratio of how many times the genealogy is updated compared to the parameters\n\
#                       If the value is 0.4 in a 2 population scenario and with 1000000 steps\n\
#                       The tree will be evaluated 400000 times, Theta_1, Theta_2, M_21, and M_12\n\
#                        will be each evaluated 125000 times.\n\
#        bayes-posteriorbins=VALUE VALUE\n\
#            VALUE      is the number of bins in the psterior distribution histogram for Theta or M\n\
#        bayes-posteriormaxtype=< ALL | P99 | MAXP99 | P100 >\n\
#            ALL        plots the WHOLE prior-parameter range\n\
\n\
#            P99        plots from the minimum prior range value to\n\
\n\
#                       the 99% percentile value of EACH parameter\n\
\n\
#            MAXP99     sets all axes from minimum to the maximal\n\
\n\
#                       99% percentile value of ALL parameter\n\
\n\
#            P100       plots from the minimum prior range value to\n\
\n\
#                       the 100% percentile value of EACH parameter\n\
\n\
#        bayes-file=<YES:FILENAME|NO>\n\
#            FILENAME is the name of the file that will contain\n\
#                    the results for the posterior distribution\n\
#        bayes-allfile=<<YES|TEMP>:INTERVAL:FILENAME|NO>\n\
#            FILENAME is the name of the file that will contain\n\
#                    all parameters of the posterior distribution [HUGE]\n\
#            INTERVAL is the interval at which all parameters are written to file\n\
\n\
#        \n\
#        bayes-proposals= THETA < SLICE | METROPOLIS >\n\
#        bayes-proposals= MIG < SLICE | METROPOLIS >\n\
#               SLICE uses the slice sampler to propose new parameter values\n\
#               METROPOLIS uses the Metropolis-Hastings sampler\n\
#               (this is done for each parameter group: THETA or MIGration)\n\
#        \n\
#        bayes-priors= THETA <UNIFORM unipriorvalues | EXP exppriorvalues | WINDOWEXP wexppriorvalues \n\
#        bayes-priors= MIG <UNIFORM unipriorvalues | EXP exppriorvalues | WINDOWEXP wexppriorvalues \n\
#                unipriorvalues: min max delta\n\
#                exppriorvalues: min mean max\n\
#                wexppriorvalues: min mean max delta\n\
#\n\
# Maximum likelihood OPTIONS\n\
#        short-chains=VALUE   VALUE is 1..n [Default is 10]\n\
#        short-inc=VALUE      VALUE is the number of updates that are not recorded\n\
#        short-sample=VALUE   VALUE is the number of sampled updates\n\
#\n\
# Search OPTIONS for both strategies\n\
#        long-chains=VALUE   VALUE is 1..n [Default is 3]\n\
#        long-inc=VALUE      VALUE is the number of updates that are not recorded\n\
#        long-sample=VALUE   VALUE is the number of sampled updates\n\
#        burn-in=VALUE       VALUE is the number of updates to discard at the beginning\n\
#        auto-tune=YES:VALUE  VALUE the the target acceptance ratio\n\
#                             if value is missing, it is set to 0.44\n\
#\n\
bayes-update=YES\n\
bayes-updatefreq=0.500000\n\
bayes-posteriorbins=1500 1500\n\
bayes-posteriormaxtype=TOTAL\n\
bayes-file=NO\n\
bayes-allfile=YES:1:bayesallfile_{0}.gz\n\
bayes-proposals= THETA METROPOLIS-HASTINGS Sampler\n\
bayes-proposals= MIG METROPOLIS-HASTINGS Sampler\n\
bayes-priors= THETA UNIFORMPRIOR: {1} {2} {3} \n\
bayes-priors= MIG UNIFORMPRIOR: {4} {5} {6} \n\
#\n\
long-chains=1\n\
long-inc={7}\n\
long-sample={8}\n\
burn-in={9}  \n\
auto-tune=YES:0.440000\n\
#\n\
'.format(run_name, str(theta_pmin), str(theta_pmax), str(theta_palpha), str(m_pmin), str(m_pmax), str(m_palpha), str(chain_inc), str(chain_sam), str(burnin)))
	parm.writelines('\
#-------------------------------------------------------------------------------\n\
# Schemes to improve MCMC searching and/or thermodynamic integration\n\
#\n\
# Heating schemes {MCMCMC = MC cubed}\n\
#    Syntax: heating=< NO | <YES | ADAPTIVE>:SKIP:TEMPERATURES\n\
#        NO    No heating\n\
#        YES   heating using TEMPERATURES\n\
#        ADAPTIVE adaptive heating using start TEMPERATURES [fails sometimes!!]\n\
#        SKIP skip that many comparisons, this lengthens the run by SKIP\n\
#            TEMPERATURES    { 1.0, temp1, temp2, temp3 .. tempn}\n\
#     Example: heating=YES:1:{1.0, 1.2, 3.0,1000000.0}\n\
# Heating:  swapping chains\n\
#     Syntax: heated-swap=< YES | NO >\n\
#         YES  swapping of chains enabled [DEFAULT]\n\
#         NO   swapping of chains disabled\n\
#      Example: heated-swap=YES\n\
heating=NO\n\
#\n\
# Lengthening chain schemes\n\
#    Syntax: moving-steps=< NO | YES:VALUE>\n\
#       VALUE   frequency is between 0..1\n\
moving-steps=NO\n\
#\n\
#    Syntax: long-chain-epsilon=VALUE\n\
#       VALUE    is between 0..INFINITY\n\
#                the VALUE is the likelihood ratio between the old and thew chain\n\
#                the VALUE depends on the number of parameters: with 1 values of 0.5 are great\n\
#                but with many parameters values and bad data >20 is more reasonable\n\
long-chain-epsilon=INFINITY\n\
#\n\
#    Convergence statistic [Gelman and Rubin]\n\
#    Syntax: gelman-convergence=< YES:Pairs|Summary | NO >\n\
#       NO      do not use Gelman\'s convergence criterium\n\
#       YES     use Gelman\'s convergence criteria between chain i, and i-1\n\
#               PAIRS reports all replicate pairs\n\
#               SUM   reports only mean and maxima\n\
gelman-convergence=No\n\
#\n\
#    Syntax: replicate=< NO | YES:<VALUE | LastChains> >\n\
#       NO     no replication of run\n\
#       YES    replicate run\n\
#           VALUE     number between 2 and many, complete replicates\n\
#           LastChains  replications over last chains\n\
replicate=NO\n\
#\n\
# Migration rates are attracted to zero (fatal attraction)\n\
# Resistance is the lowest migration value for all but the last chain\n\
#    Syntax resistance=VALUE\n\
#        VALUE is the lowest migration rate value allowed during all but the last chain\n\
#              typical values are 0.01 or _lower_ for data with sequences and 0.0001 or _lower_ for other data\n\
resistance=0.000001\n\
#\n\
#-------------------------------------------------------------------------------\n\
#\n\
end\n\
')
	
	parm.close() # close the parmfile