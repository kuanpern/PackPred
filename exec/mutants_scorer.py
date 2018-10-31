import os
import sys
import json
import argparse
from numpy import *
import cPickle as pickle
from copy import deepcopy
from modeller import *
log.none()
env = environ()

from clique_scoring_helper import *

# read inputs
parser = argparse.ArgumentParser(description='Score cliques')
parser.add_argument('--clique', type=str, required=True, help="input clique file")
parser.add_argument('--pdb', type=str, required=True, help="input pdb file")
parser.add_argument('--mutation', type=str, required=True, help="mutation list file")
parser.add_argument('--scoredir', type=str, required=True, help="scores source directory")
parser.add_argument('--outfile', type=str, required=True, help="output file name")
parser.add_argument('--dcut',  type=float, required=False, help="distance cut-off (Angstrom)", default=10.0)
args = parser.parse_args()
pars = vars(args)

# parse inputs
infilename          = pars['clique']
pdb_fname           = pars['pdb']
mutation_list_fname = pars['mutation']
score_files_dir     = pars['scoredir']
outfname            = pars['outfile']
dcut                = pars['dcut']

d_cuts = [dcut]
d_cuts.sort()

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


