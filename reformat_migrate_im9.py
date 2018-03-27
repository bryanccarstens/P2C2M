##### Script to use seq-gen data to create Migrate (datafile and paramter file), IMa, and fasta files #####

### Import libraries
import glob # import glob module used to obtain file names
from natsort import natsorted # Note: must install natsort module for python 3 "sudo pip3 install natsort". This performs a natural sort, i.e. (1,2,3,4,5,6,7,8,9,10) rather than (1,10,2,3,4,5,6,7,8,9)
import re # import regular expressions library for creating fasta files
from write_parm import write_parm # imports function to write the migrate parameter files. Note: requires write_parm.py script in same directory. Migrate settings can be adjusted in this file

### Define values
num_pops = 3 # number of populations. Note: changing the number of pops here does not change everything. See lines 18-21, 61-62
num_inds = 20 # number alleles/individuals per population
num_loci = 10 # number of loci
popind = (num_pops * num_inds) +1 # number of lines making up each locus in the seq-gen file. used for counter below

### Get information from seq-gen file and write new IMa and migrate files
seq_names = glob.glob('seqfile*') # obtain names for all seq-gen output files

for file_name in seq_names: # for each seq-gen file
	name_cat = file_name.split('_')[1] + '_' + file_name.split('_')[2] + '_' + file_name.split('_')[3] # remove seqfile from name
	
	#Note: Add more lines here when adding populations
	popA = open('popA', 'w+') # create popA file 
	popB = open('popB', 'w+') # create popB file
	popC = open('popC', 'w+') # create popC file
	
	seqgen = open(file_name, 'rt') # open seq-gen file
	seq_file = seqgen.readlines() # read all lines of seq-gen file, each line becomes item in list seq_file
	seqgen.close() # close connection to seq-gen file
	
	### IMa and fasta
	ima = open('ima_' + name_cat, 'at') # create IMa file
	
	# Note: if num_pops is adjusted, the pop names and tree will need to be adjusted in line below
	ima.writelines(name_cat + '\n' + str(num_pops) + '\n' + 'popA popB popC\n' + '((0,1):3,2):4\n' + str(num_loci) +'\n') # write first 5 lines to IMa file 
	
	s = 0 # counter start. Always start at 0
	e = popind # counter end. number of lines making up each locus in the seq-gen file
	for locus in range(1, num_loci+1): # for each locus
		loc_lines = seq_file[s:e] # subset seq-gen lines for locus 
		loc_lines = loc_lines[1:] # remove locus info line
		loc_lines = natsorted(loc_lines) # sort lines by pop/ind
			
		### IMa
		ima.writelines('locus' + str(locus) + ' ' + ((str(num_inds) + ' ')*num_pops) + '1000 H 1.0\n') # write first line for locus in IMa file
		
		ima.writelines(loc_lines) # write sorted lines to IMa file
		
		# Note: add more lines here if increasing num_pops
		popA.writelines(loc_lines[0:num_inds]) # write popA data to popA file. 
		popB.writelines(loc_lines[num_inds:(2*num_inds)]) # write popB data to popB file
		popC.writelines(loc_lines[(2*num_inds):]) # write popC data to popC file
		
		### fasta
		fasta_name = '{0}_locus{1}.fasta'.format(name_cat, locus) # create fasta file name for locus
		fasta = open(fasta_name, 'w') # create fasta file
		for sample in loc_lines: # for each individual
			line = re.sub(r'(pop)', r'>\1', sample) # add '>' for sample names
			line = re.sub(r'\s+', r'\n', line) # move sequence to new line below sample name
			line1 = line.split('\n')[0] # separate sample name
			line2 = line.split('\n')[1] # separate sequence
			line2 = '\n'.join(line2[i:i+80] for i in range(0, len(line2), 80)) # create new line for every 80bp
			fasta.writelines(line1 + '\n' + line2 + '\n') # write fasta entry to fasta file
		fasta.close() # close fasta file
		
		s+=popind # adjust counter for next locus
		e+=popind # adjust counter for next locus
	
	ima.writelines('\n') # add empty line to end of IMa file
	ima.close() # close IMa file

	### migrate-n
	write_parm(name_cat) # write the migrate parameter file
	
	migrate = open('migrate_' + name_cat, 'at') # create migrate file
	migrate.writelines(str(num_pops) + ' ' + str(num_loci) + ' ' + name_cat + '\n' + ('1000 '*num_loci) + '\n') # write first 2 lines to migrate file
	
	# Note: add more names to line below if increasing num_pops
	pop_names = ['popA', 'popB', 'popC'] # create list of population names. 
	for pop in pop_names:
		migrate.writelines(((str(num_inds) + ' ') * num_loci) + pop + '\n') # write population identifier line
		eval(pop).seek(0) # move cursor back to beginning of pop file so lines can be read in
		pop_lines = eval(pop).readlines() # read all lines of population file, each line becomes item in list pop_lines
		
		b = 0 # counter start. Always start at 0
		t = num_inds # counter end. number of lines making up each locus in the pop text files
		for locus in range(1, num_loci+1): # for each locus
			loc_lines = pop_lines[b:t] # subset pop lines for locus
			migrate.writelines(loc_lines) # write sorted lines to migrate file
			migrate.writelines('\n') # add empty line 
			
			b+=num_inds # adjust counter for next locus
			t+=num_inds # adjust counter for next locus
			
	migrate.close() # close migrate file
		
		 