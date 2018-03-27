### Script for running migrate-n over many files ###
# Note: this script assumes all parmfiles are in the same directory as the migrate-n executable and all parmfiles begin with 'parmfile'

# Import Libraries
import glob # import the glob module to obtain file names
import subprocess # import the subprocess module for running other programs

### Get parameter file names
parm_names = glob.glob('parmfile*') # obtain names for all migrate parameter files
for file_name in parm_names: # for each parameter file
	cmd_line = './migrate-n {0} -nomenu'.format(file_name) # create string for command to run migrate-n
	subprocess.call(cmd_line, shell=True) # run migrate-n
	name_cat = file_name.split('_')[1] + '_' + file_name.split('_')[2] + '_' + file_name.split('_')[3] # remove parmfile from name
	cmd_line2 = 'mv treefile migrate_treefile_{0}'.format(name_cat) # create string for command to rename treefile
	subprocess.call(cmd_line2, shell=True) # change name of treefile