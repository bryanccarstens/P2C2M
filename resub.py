##### Function for replacing numbers in ms taxon identifiers with population name and individual number #####
# This script is currently used with the sim_migrate_im3.py script for simulating data to use in migrate-n and IMa
# file_name refers to a file name containing the ms gene trees
# If the number of individuals per population changes, adjust the numbers in lines 16-20
#	Numbers in brackets indicate a range for the ones digit. The number preceding the brackets is the tens digit
#	Currently, we are using 20 inds/pop, so popA is inds 1-20, popB is inds 21-40, and popC is inds 41-60
# If the number of simulated populations changes, add lines similar to lines 16-20, changing the 'pop' name
#	Also adjust the numbers to match the populations as described above

import re # import the regular expressions library

def tax_names(file_name): # create tax_names function that requires input of file_name
	tree = open(file_name, 'r') # open ms file just created
	tree_lines = tree.readlines() # read in lines from ms treefile. Each line is a separate gene tree
	tree.close() # close ms treefile
	treefile = open(file_name, 'w+') # overwrite the ms treefile
	for line in tree_lines: # for each line/tree
		line = re.sub(r'\(([1-9]:|1[0-9]:|2[0]:)', r'(popA\1', line) # change taxon identifiers 1-20 to popA1-popA20
		line = re.sub(r',([1-9]:|1[0-9]:|2[0]:)', r',popA\1', line) # change taxon identifiers 1-20 to popA1-popA20
		line = re.sub(r'\((2[1-9]:|3[0-9]:|4[0]:)', r'(popB\1', line) # change taxon identifiers 21-40 to popB21-popB40
		line = re.sub(r',(2[1-9]:|3[0-9]:|4[0]:)', r',popB\1', line) # change taxon identifiers 21-40 to popB21-popB40
		line = re.sub(r'\((4[1-9]:|5[0-9]:|6[0]:)', r'(popC\1', line) # change taxon identifiers 41-60 to popC41-popC60
		line = re.sub(r',(4[1-9]:|5[0-9]:|6[0]:)', r',popC\1', line) # change taxon identifiers 41-60 to popC41-popC60
		treefile.writelines(line) # write corrected lines to treefile
	
	treefile.close() # close the new treefile
	
