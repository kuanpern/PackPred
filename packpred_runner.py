#!/usr/bin/python

# 1. generate mutations
# 2. compute depth
# 3. compute clique
# 4. perform mutation
# 5. score mutation
# 6. collect scores; charting, if necessary / requested

import os
import sys
import time
import json
import uuid
import glob
import commands
import argparse
import cPickle as pickle
from string import Template
from scipy.interpolate import *

# setup logging facilities
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# 0. read input
parser = argparse.ArgumentParser(description='Packpred program')
parser.add_argument('--input', type=str, required=True, help="input parameters (.json)")
parser.add_argument('--dcut',  type=float, required=False, help="distance cut-off (Angstrom)", default=10.0)
parser.add_argument('--order', type=int, required=False, help="maximum clique order", choices=(2, 3, 4, 5), default=3)
parser.add_argument('--exedir', type=str, required=True, help="executable directory")
parser.add_argument('--workdir', type=str, required=True, help="working directory")
pars = parser.parse_args()
pars = vars(parser)

# parse inputs
d_cut            = pars['dcut']
max_clique_order = pars['order']
input_json_fname = pars['input']
exedir           = pars['exedir']
workdir          = pars['workdir']

# set parameters
min_clique_order = 2
max_clique_cutoff = max(d_cuts) + 0.5
d_cuts = [d_cut]

# read json input
with open(input_json_fname) as fin:
	inputs = json.loads(fin.read())
# end with

# retrieve job-id for file naming
if 'job_id' in inputs.keys():
	job_id = str(inputs['job id'])
else:
	job_id = str(uuid.uuid4())
# end if

# setup the logging facility
if 'log filename' in inputs:
	logfilename = inputs['log filename']
	# create a file handler
	handler = logging.FileHandler(logfilename)
else:
	handler = logging.StreamHandler()
# end if
# create a logging format
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# writing to log and err files
def stream_to_log(msg):
	assert type(msg) in [list, tuple, str, unicode], 'message type "'+str(type(msg))+'" not understood'
	if type(msg) in [list, tuple]:
		msg = [str(t) for t in msg]
		msg = ' '.join(msg) + '\n'
	# end if

	logger.info(msg)
# end def

def stream_to_err(msg, exit_code = 1):
	assert type(msg) in [list, tuple, str, unicode], 'message type "'+str(type(msg))+'" not understood'
	if type(msg) in [list, tuple]:
		msg = [str(t) for t in msg]
		msg = ' '.join(msg) + '\n'
	# end if

	logger.error(msg)
	# exit the program also
	sys.exit(exit_code)
# end def

# -- error checking - (1)
if not(os.path.isfile('native.pdb')):
	stream_to_err('expecting an "native.pdb" file, not found')
# end if

stream_to_log(' - reading input ...' )
# parse input
# input_pdb = str(inputs['pdb name'])
# pdb_code  = input_pdb.split('.pdb')[0]
input_pdb = 'native.pdb'
pdb_code  = 'native'
	
# 1. generate mutation list
stream_to_log(' - generate mutation list')
# 1.0. specific program path
mutation_list_fname = 'mutation.list'
generate_mutation_exec = exedir+'/exec/generate_mutations.py'
cmd = ['python', generate_mutation_exec, input_json_fname, mutation_list_fname]
cmd = ' '.join(cmd)
stream_to_log([ cmd, '\n', commands.getoutput(cmd) ])

# -- error checking - (2)
if not(os.path.isfile(mutation_list_fname)):
	stream_to_err('mutation list file: '+mutation_list_fname+' not found')
# end if
if len(open(mutation_list_fname).read().splitlines()) == 0:
	stream_to_err('mutation file empty')
# end if

# 2. compute depth
stream_to_log(' - computing residue depths ...')
# 2.1. invoke depth program
depth_exec = exedir+'/exec/depth-2.0/bin/DEPTH'
cmd = [depth_exec, '-i', input_pdb, '-o', input_pdb, '-survive 4']
cmd = ' '.join(cmd)
stream_to_log([ cmd, '\n', commands.getoutput(cmd) ])
# 2.2. move files into depth directory
os.makedirs('depths/')
cmd = 'mv native*depth* depths/'
stream_to_log([ cmd, '\n', commands.getoutput(cmd) ])

# -- error checking - (3)
depth_fname = 'depths/'+input_pdb+'-residue.depth'
if not(os.path.isfile(depth_fname)):
	stream_to_err('depth file: '+depth_fname+' not found')
# end if

# 3. compute clique
stream_to_log(' - Computing clique ...')
# 3.0. specify clique program paths
clique_exedir	         = exedir + '/exec/YETI-2.0/'
clique_main_exec		 = clique_exedir + '/main.py'
clique_distance_selector = clique_exedir + '/exec/distance_selector.py'
clique_linkage_builder   = clique_exedir + '/exec/linkage'
clique_linkage_par_file  = clique_exedir + '/lib/half_side.density'
clique_identifier_exe	 = clique_exedir + '/exec/identify_clique.py'

# 3.1. invoke clique program
cmd = [exedir+'/venv/bin/python ', clique_main_exec, input_pdb, str(max_clique_order), str(max_clique_cutoff), clique_distance_selector, clique_linkage_builder, clique_linkage_par_file]
cmd = ' '.join(cmd)
stream_to_log([cmd, '\n', commands.getoutput(cmd) ])

# 3.2. move files into respective directories
cmds = Template('''
mkdir cliques/ tempfiles/
mv $pdb_code*.clique? cliques/
mv $pdb_code*.link $pdb_code*orphan* $pdb_code*xyz $pdb_code*select.pdb $pdb_code*distmap tempfiles/
''').substitute({'pdb_code': pdb_code})
cmds = cmds.strip().splitlines()
for cmd in cmds:
	stream_to_log([cmd, '\n', commands.getoutput(cmd) ])
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

	cmd = [exedir+'/venv/bin/python', clique_identifier_exe, clique_filename, pdb_fname, depth_fname, output_filename]
	cmd = ' '.join(cmd)
	stream_to_log([cmd, '\n', commands.getoutput(cmd) ])
# end for

# -- error checking - (4)
clique_fname = 'native.clique3.pickle'
if not(os.path.isfile(clique_fname)):
	stream_to_err('essential clique file: '+clique_fname+' not found')
# end if

# 5. score native clique
stream_to_log(' - Scoring native cliques ...')
scores_lib_dir	    = exedir + '/lib/potentials/'
clique_scorer_exec  = exedir + '/exec/clique_scorer.py'
# 5.1. scoring all cliques
clique_files = glob.glob('*.clique*.pickle')
for clique_filename_pickle in clique_files:
	cmd = [exedir+'/venv/bin/python', clique_scorer_exec, clique_filename_pickle, scores_lib_dir, '../job.specifications']
	cmd = ' '.join(cmd)
	stream_to_log([cmd, '\n', commands.getoutput(cmd) ])
# end for

# 5.2. score all mutants
stream_to_log(' - Scoring mutants ...')
mutants_scorer_exec  = exedir + '/exec/mutants_scorer.py'
clique_score_fnames = []
for clique_filename_pickle in clique_files:
	clique_score_outfname = clique_filename_pickle + '.json'
	clique_score_fnames.append(clique_score_outfname)

	cmd = [exedir+'/venv/bin/python', mutants_scorer_exec, clique_filename_pickle, '../'+input_pdb, '../'+mutation_list_fname, scores_lib_dir, clique_score_outfname, '../job.specifications']
	cmd = ' '.join(cmd)
	stream_to_log(cmd)
	output = commands.getoutput(cmd).splitlines()
# end for

# record scores from different orders
Scores = {}
for clique_score_outfname in clique_score_fnames:
	with open(clique_score_outfname) as fin:
		mutant_score_data = json.load(fin)
		Scores.update(mutant_score_data)
	# end with
# end for

# 5.3. combination formula - just single score ...
combined_scores = {}
with open('../'+mutation_list_fname) as fin:
	mutation_keys = json.load(fin)['mutations']
# end with
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
	stream_to_err('no mutation scorings were performed.')
# end if

# 6. generate html output
# 6.1. tsv output
os.chdir(cwd)
tmpfname = job_id + '.tmp.json'
with open(tmpfname, 'w') as fout:
	json.dump(fout, combined_scores)
# end with

# 6.2. generate html output
stream_to_log(' - generating output html')
result_generator_exec = exedir + '/exec/result_generator.py'
output_html_fname = inputs['output filename']
template_fname = exedir + '/templates/output_template.html'
colorbar_fname = exedir + '/templates/colorbar64.txt'
cmd = [exedir+'/venv/bin/python', result_generator_exec, tmpfname, template_fname, colorbar_fname, output_html_fname]
cmd = ' '.join(cmd)
stream_to_log([cmd, '\n', commands.getoutput(cmd) ])


