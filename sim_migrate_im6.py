##### Note on abbreviations and values #####
# When possible, abbreviations follow those defined in the ms manual. M refers to 4Nom, where No is the effective population size and m is the fraction of each subpopulation made up of new migrants each generation. 
# Divergence times t1 and t2 are in units of 4No generations. 


### Import Libraries
import random # import the random library
import subprocess # import the subprocess module for running other programs
import time # import the timing module for measuring how long it takes the script to run
from resub import tax_names # import function for changing taxonomy names in the ms treefiles. Note: Requires resub.py script in same directory. Items in this function will need to be adjusted in order to change the number of populations/inds per population.

t = time.process_time() # start timer

### Define Values
reps = 10 # number of reps to perform
M_a = 100 # parameter value used to specify migration distribution. See the different distributions below for determining a value of a
M_b = 10000 # parameter value used to specify migration distribution. See the different distributions below for determining a value of b
t1_a = 1.5 # parameter value used to specify divergence time 1 distribution. See the different distributions below for determining a value of a. Units = 4No generations
t1_b = 2.5 # parameter value used to specify divergence time 1 distribution. See the different distributions below for determining a value of b. Units = 4No generations
t2_a = 0.5 # parameter value used to specify divergence time 2 distribution. See the different distributions below for determining a value of a. Units = 4No generations
t2_b = 1.5 # parameter value used to specify divergence time 2 distribution. See the different distributions below for determining a value of b. Units = 4No generations

sim_log = open('sim_log.txt', 'w') # create a log file to record the migration and divergence time values chosen from the distributions for each simulation rep
sim_log.writelines('Rep,Mig,t1,t2,t2_sec') # write header line
##### Distributions #####
# For each distribution, there are multiple types to choose from. In the lines below, they are specified by 'random.x(a, b)' where x is the distribution type.
# 	To change the distribution used, adjust lines 45-47 accordingly:
#		uniform - a uniform distribution. a and b correspond to the min and max
#		betavariate - a beta distribution. a and b correspond to the alpha and beta parameters
#		expovariate - an exponential distribution. This distribution will take the 'a' value only (i.e. random.expovariate(a))! a corresponds to 1.0 divided by the desired mean
#		gammavariate - a gamma distribution. a and b correspond to the alpha and beta parameters
#		gauss - a Gaussian distribution. a and b correspond to the mean and standard deviation
#		normalvariate - a normal distribution. a and b correspond to the mean and standard deviation

##### ms and seq-gen input values #####
# Currently ms is run with the values under the define values section of this script along with the following:
#	number of populations = 3
#	number of individuals per population = 20
#	number of loci = 10
#	theta = 1 per locus (divided by 1000bp per locus = 0.0001)
#	mutation model (for seq-gen) = HKY
#	branch length scaling factor = 0.1
# See ms and seq-gen manuals for more information

### Conduct Simulations
for i in range(1, reps + 1): # over how many reps. Note: +1 is used because range does not include the last number
	M = random.uniform(M_a, M_b) # select a migration rate at random from a uniform distribution. Units = 4Nom
	t1 = random.uniform(t1_a, t1_b) # select a divergence time at random from a uniform distribution. Units = 4No generations
	t2 = random.uniform(t2_a, t2_b) # select a divergence time at random from a uniform distribution. Units = 4No generations
	t2_secondary = 0.5*t2 # calculate time of secondary contact. Units = 4No generations
	
	sim_log.writelines(str(i) + ',' + str(M) + ',' + str(t1) + ',' + str(t2) + ',' + str(t2_secondary) + '\n') # write parameter values to log
	
	### Run ms and seq-gen for each model
	# n-island model
	msmod1 = './ms 60 10 -t 1.0 -I 3 20 20 20 {0} -T | tail -n +4 | grep -v //>treefile_model1_rep{1}'.format(M, i) # create string for command to run ms
	subprocess.call(msmod1, shell=True) # run ms
	file_name1 = 'treefile_model1_rep{}'.format(i) # create file name to be used with tax_names function
	tax_names(file_name1) # run function to change taxon names
	sgmod1_20 = './seq-gen -mHKY -l 1000 -s 0.0016 <treefile_model1_rep{0} > seqfile_model1_rep{0}_s20'.format(i) # create string for command to run seq-gen with 20 segregating sites
	sgmod1_60 = './seq-gen -mHKY -l 1000 -s 0.0047 <treefile_model1_rep{0} > seqfile_model1_rep{0}_s60'.format(i) # create string for command to run seq-gen with 60 segregating sites
	subprocess.call(sgmod1_20, shell=True) # run seq-gen
	subprocess.call(sgmod1_60, shell=True) # run seq-gen
	
	# continuous migration
	msmod2 = './ms 60 10 -t 1.0 -I 3 20 20 20 {0} -ej {1} 2 1 -ej {2} 1 3 -T | tail -n +4 | grep -v //>treefile_model2_rep{3}'.format(M, t2, t1, i) # create string for command to run ms
	subprocess.call(msmod2, shell=True) # run ms
	file_name2 = 'treefile_model2_rep{}'.format(i) # create file name to be used with tax_names function
	tax_names(file_name2) # run function to change taxon names
	sgmod2_20 = './seq-gen -mHKY -l 1000 -s 0.0016 <treefile_model2_rep{0} > seqfile_model2_rep{0}_s20'.format(i) # create string for command to run seq-gen with 20 segregating sites
	sgmod2_60 = './seq-gen -mHKY -l 1000 -s 0.0047 <treefile_model2_rep{0} > seqfile_model2_rep{0}_s60'.format(i) # create string for command to run seq-gen with 60 segregating sites
	subprocess.call(sgmod2_20, shell=True) # run seq-gen
	subprocess.call(sgmod2_60, shell=True) # run seq-gen
	
	# secondary contact
	msmod3 = './ms 60 10 -t 1.0 -I 3 20 20 20 {0} -ej {1} 2 1 -ej {2} 1 3 -em {3} 2 1 {0} -T | tail -n +4 | grep -v //>treefile_model3_rep{4}'.format(M, t2, t1, t2_secondary, i) # create string for command to run ms
	subprocess.call(msmod3, shell=True) # run ms
	file_name3 = 'treefile_model3_rep{}'.format(i) # create file name to be used with tax_names function
	tax_names(file_name3) # run function to change taxon names
	sgmod3_20 = './seq-gen -mHKY -l 1000 -s 0.0016 <treefile_model3_rep{0} > seqfile_model3_rep{0}_s20'.format(i) # create string for command to run seq-gen with 20 segregating sites
	sgmod3_60 = './seq-gen -mHKY -l 1000 -s 0.0047 <treefile_model3_rep{0} > seqfile_model3_rep{0}_s60'.format(i) # create string for command to run seq-gen with 60 segregating sites
	subprocess.call(sgmod3_20, shell=True) # run seq-gen
	subprocess.call(sgmod3_60, shell=True) # run seq-gen
	
	# isolation only
	msmod4 = './ms 60 10 -t 1.0 -I 3 20 20 20 0.0 -ej {0} 2 1 -ej {1} 1 3 -T | tail -n +4 | grep -v //>treefile_model4_rep{2}'.format(t2, t1, i) # create string for command to run ms
	subprocess.call(msmod4, shell=True) # run ms
	file_name4 = 'treefile_model4_rep{}'.format(i) # create file name to be used with tax_names function
	tax_names(file_name4) # run function to change taxon names
	sgmod4_20 = './seq-gen -mHKY -l 1000 -s 0.0016 <treefile_model4_rep{0} > seqfile_model4_rep{0}_s20'.format(i) # create string for command to run seq-gen with 20 segregating sites
	sgmod4_60 = './seq-gen -mHKY -l 1000 -s 0.0047 <treefile_model4_rep{0} > seqfile_model4_rep{0}_s60'.format(i) # create string for command to run seq-gen with 60 segregating sites
	subprocess.call(sgmod4_20, shell=True) # run seq-gen
	subprocess.call(sgmod4_60, shell=True) # run seq-gen

sim_log.close() # close the log file	
elapsed_t = time.process_time() - t # calculate time elapsed
print('\nTime elapsed: ', elapsed_t, ' seconds') # print time elapsed
	
