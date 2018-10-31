#!/usr/bin/python
# 'one-stop' script to compute clique score
# input: input parameters in pickle format

# 1. generate mutations
# 2. compute depth
# 3. compute clique
# 4. perform mutation
# 5. score mutation
# 6. collect scores; charting, if necessary / requested

import sys
import time
import commands
import os
import json
import pickle
from string import Template
from scipy.interpolate import *

# - install directory and paths
install_dir = '/home/tankp/Desktop/Mutational_stability_prediction/Content/programs/yeti_mutational_perdictor_1.2/'

# = computation starts here = #

# 0. read input
if len(sys.argv) != 2:
	print 'python', sys.argv[0], 'input.json'
	sys.exit(0)
# end if

input_json_fname = sys.argv[1]
inputs = json.loads(open(input_json_fname).read())

# retrieve job-id for file naming
job_id = str(inputs['job id'])

# --- check code --- #
logfilename = inputs['log filename']
errfilename = inputs['err filename']

cwd = os.getcwd()
errfilename = cwd + '/' + os.path.basename(errfilename)
logfilename = cwd + '/' + os.path.basename(logfilename)

# writing to log and err files
def print_to_log(filename, msg):
	assert type(msg) in [list, tuple, str, unicode]
	if type(msg) in [list, tuple]:
		msg = [str(t) for t in msg]
		msg = ' '.join(msg) + '\n'

	open(filename, 'a').writelines(str(msg))
# end def

def print_to_error(filename, msg, exit_code = 1):
	open(filename, 'a').writelines(msg)
	sys.exit(exit_code)
# end def

# -- error checking - (1)
if not(os.path.isfile('native.pdb')):
	print_to_error(errfilename, 'expecting an "native.pdb" file, not found')
# end if

# copy necessary files into the directory
cmd = 'ln -sf '+install_dir+'/lib/specifications job.specifications'
print_to_log(logfilename, [ cmd, '\n', commands.getoutput(cmd) ])
# read specifications
execfile('job.specifications') # <- exec(), you sure ?

print_to_log(logfilename, [ ' - reading input ...' ])
# parse input
# input_pdb = str(inputs['pdb name'])
# pdb_code  = input_pdb.split('.pdb')[0]
input_pdb = 'native.pdb'
pdb_code  = 'native'
	
# 1. generate mutation list
print_to_log(logfilename, [ ' - generate mutation list' ])
# 1.0. specific program path
mutation_list_fname = 'mutation.list'
generate_mutation_exec = install_dir+'/exec/generate_mutations.py'
cmd = ['python', generate_mutation_exec, input_json_fname, mutation_list_fname]
cmd = ' '.join(cmd)
print_to_log(logfilename, [ cmd, '\n', commands.getoutput(cmd) ])

# -- error checking - (2)
if not(os.path.isfile(mutation_list_fname)):
	print_to_error(errfilename, 'mutation list file: '+mutation_list_fname+' not found')
if len(open(mutation_list_fname).read().splitlines()) == 0:
	print_to_error(errfilename, 'mutation file empty')

# 2. compute depth
print_to_log(logfilename, [ ' - computing residue depths ...' ])
# 2.1. invoke depth program
depth_exec = install_dir+'/exec/depth-2.0/bin/DEPTH'
cmd = [depth_exec, '-i', input_pdb, '-o', input_pdb, '-survive 4']
cmd = ' '.join(cmd)
print_to_log(logfilename, [ cmd, '\n', commands.getoutput(cmd) ])
# 2.2. move files into depth directory
os.makedirs('depths/')
cmd = 'mv native*depth* depths/'
print_to_log(logfilename, [ cmd, '\n', commands.getoutput(cmd) ])

# -- error checking - (3)
depth_fname = 'depths/'+input_pdb+'-residue.depth'
if not(os.path.isfile(depth_fname)):
	print_to_error(errfilename, 'depth file: '+depth_fname+' not found')

# 3. compute clique
print_to_log(logfilename, [' - Computing clique ...' ])
# 3.0. specify clique program paths
clique_install_dir	   = install_dir + '/exec/YETI-2.0/'
clique_main_exec		 = clique_install_dir + '/main.py'
clique_distance_selector = clique_install_dir + '/exec/distance_selector.py'
clique_linkage_builder   = clique_install_dir + '/exec/linkage'
clique_linkage_par_file  = clique_install_dir + '/lib/half_side.density'
clique_identifier_exe	= clique_install_dir + '/exec/identify_clique.py'

# 3.1. invoke clique program
cmd = ['python ', clique_main_exec, input_pdb, str(max_clique_order), str(max_clique_cutoff), clique_distance_selector, clique_linkage_builder, clique_linkage_par_file]
cmd = ' '.join(cmd)
print_to_log(logfilename, [cmd, '\n', commands.getoutput(cmd) ])

# 3.2. move files into respective directories
cmds = Template('''
mkdir cliques/ tempfiles/
mv $pdb_code*.clique? cliques/
mv $pdb_code*.link $pdb_code*orphan* $pdb_code*xyz $pdb_code*select.pdb $pdb_code*distmap tempfiles/
''').substitute({'pdb_code': pdb_code})
cmds = cmds.strip().splitlines()
for cmd in cmds:
	print_to_log(logfilename, [cmd, '\n', commands.getoutput(cmd) ])
# end for

# 3.3. identify residue clique (append residue / environment information to clique file)
cwd = os.getcwd()
os.chdir('cliques')
clique_files = commands.getoutput('ls *.clique*').splitlines()
for clique_filename in clique_files:
	root = os.path.basename(clique_filename).split('.clique')[0]
	pdb_fname   = '../'+root+'.pdb'
	depth_fname = '../depths/'+input_pdb+'-residue.depth'
	output_filename = clique_filename+'.pickle'

	cmd = ['python', clique_identifier_exe, clique_filename, pdb_fname, depth_fname, output_filename]
	cmd = ' '.join(cmd)
	print_to_log(logfilename, [cmd, '\n', commands.getoutput(cmd) ])
# end for

# -- error checking - (4)
clique_fname = 'native.clique3.pickle'
if not(os.path.isfile(clique_fname)):
	print_to_error(errfilename, 'essential clique file: '+clique_fname+' not found')

# 5. score native clique
print_to_log(logfilename, [' - Scoring native cliques ...' ])
scores_lib_dir	  = install_dir + '/lib/potentials/'
clique_scorer_exec  = install_dir + '/exec/clique_scorer.py'
# 5.1. scoring all cliques
clique_files = commands.getoutput('ls *.clique*.pickle').splitlines()
for clique_filename_pickle in clique_files:
	cmd = ['python', clique_scorer_exec, clique_filename_pickle, scores_lib_dir, '../job.specifications']
	cmd = ' '.join(cmd)
	print_to_log(logfilename, [cmd, '\n', commands.getoutput(cmd) ])
# end for

# 5.2. score all mutants
print_to_log(logfilename, [' - Scoring mutants ...' ])
mutants_scorer_exec  = install_dir + '/exec/mutants_scorer.py'
clique_score_fnames = []
for clique_filename_pickle in clique_files:
	clique_score_outfname = clique_filename_pickle + '.json'
	clique_score_fnames.append(clique_score_outfname)

	cmd = ['python', mutants_scorer_exec, clique_filename_pickle, '../'+input_pdb, '../'+mutation_list_fname, scores_lib_dir, clique_score_outfname, '../job.specifications']
	cmd = ' '.join(cmd)
	print_to_log(logfilename, [cmd ])
	output = commands.getoutput(cmd).splitlines()
# end for

# record scores from different orders
Scores = {}
for clique_score_outfname in clique_score_fnames:
	mutant_score_data = json.loads(open(clique_score_outfname).read())
	Scores.update(mutant_score_data)
# end for

# 5.3. combination formula - just single score ...
combined_scores = {}
mutation_keys = json.loads(open('../'+mutation_list_fname).read())['mutations']
mutation_keys = sorted(list(set(mutation_keys) - set(['native'])))
for mut_key in mutation_keys:
	chain, mut = mut_key.split(':')
	resnum = mut[1:-1]
	resname = chain + resnum

	score = float(Scores['3'][mut_key]['10.0']) - float(Scores['3']['native']['10.0'])

	combined_scores[mut_key] = score
# end for

# -- error checking - (5)
if len(combined_scores) == 0:
	print_to_error(errfilename, 'no mutation scorings were performed.')

# 6. generate html output
# 6.1. tsv output
os.chdir(cwd)
tmpfname = job_id + '.tmp.json'
open(tmpfname, 'w').writelines(json.dumps(combined_scores))

# 6.2. generate html output
print_to_log(logfilename, [' - generating output html' ])
result_generator_exec = install_dir + '/exec/result_generator.py'
output_html_fname = inputs['output filename']
template_fname = install_dir + '/templates/output_template.html'
colorbar_fname = install_dir + '/templates/colorbar64.txt'
cmd = ['python', result_generator_exec, tmpfname, template_fname, colorbar_fname, output_html_fname]
cmd = ' '.join(cmd)
print_to_log(logfilename, [cmd, '\n', commands.getoutput(cmd) ])


