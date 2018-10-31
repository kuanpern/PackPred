import os
import sys
import json
from numpy import *
import cPickle as pickle
from copy import deepcopy
from modeller import *
log.none()
env = environ()

from clique_scoring_helper import *

# read input
if len(sys.argv) != 7:
	print 'python', sys.argv[0], 'input.cliqueN.pickle input.pdb mutation_list_file score_files_dir output.json specification_file'
	sys.exit(0)
	
infilename          = sys.argv[1]
pdb_fname           = sys.argv[2]
mutation_list_fname = sys.argv[3]
score_files_dir     = sys.argv[4]
outfname            = sys.argv[5]
specification_file  = sys.argv[6]

# import parameters
execfile(specification_file)

# read pdb file
mdl = model(env, file = pdb_fname)

# read cliques
Data = pickle.load(open(infilename))
Residues, Cliques = [Data['Residues'], Data['Cliques']]
order = Cliques[0]['order']

# - scoring of native structure
for res in Residues.values():
	score_residue(res, Cliques, d_cuts)
# end for
native_scores = score_structure(Residues, d_cuts)
Output = {'native': native_scores}

# read mutations
mutations = json.loads(open(mutation_list_fname).read())['mutations']

# read scores library
Scores_lib = read_Scores_lib(score_files_dir, d_cuts, order)

# score for every mutation
import time
start_time = time.time()
for mutation in mutations:
	Residues_mut = deepcopy(Residues)

	# get mutation information
	chain, mut = mutation.split(':')
	_from, _num, _to = mut[0], mut[1:-1], mut[-1]
	resname = chain+_num
	
	# get the mutated residue
	res = Residues_mut[resname]

	# go to all its parent cliques
	clique_indices = res['clique_memberships']
	cliques_native = {}
	for index in clique_indices:
		clique = Cliques[index]
		cliques_native[index] = deepcopy(clique) # record to put back later

		# mutate the clique and get the score
		new_clique = deepcopy(clique)
		new_clique['clique type'] = clique['clique type'].replace(_from, _to, 1)
		new_clique['clique type'] = ''.join(list(sorted(list(new_clique['clique type']))))
		new_clique = score_clique(new_clique, d_cuts, Scores_lib)
		Cliques[index] = new_clique
	# end for
	
	# update score for only residues involved in mutated cliques
	# 1. find all residues involved
	mutated_residues_index = []
	for index in clique_indices:
		mutated_residues_index += Cliques[index]['residues']
	# end for
	mutated_residues_index = list(set(mutated_residues_index))

	# 2. update their scores
	for index in mutated_residues_index:
		score_residue(Residues_mut[index], Cliques, d_cuts)
	# end for

	Output[mutation] = score_structure(Residues_mut, d_cuts)

	# revert to native Cliques
	for index in clique_indices:
		Cliques[index] = cliques_native[index]
	# end for		
# end for
end_time = time.time()

Output = {order: Output}

# write output
open(outfname, 'w').writelines(json.dumps(Output))


