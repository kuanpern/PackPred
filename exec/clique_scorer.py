import cPickle as pickle
import sys
import json
import argparse
from scipy.stats import norm
from copy import deepcopy
from numpy import *
from clique_scoring_helper import *

# read input
parser = argparse.ArgumentParser(description='Score cliques')
parser.add_argument('--infile', type=str, required=True, help="input clique file")
parser.add_argument('--scoredir', type=str, required=True, help="scores source directory")
parser.add_argument('--dcut',  type=float, required=False, help="distance cut-off (Angstrom)", default=10.0)
pars = parser.parse_args()
pars = vars(parser)

# parse inputs
clique_fname     = pars['infile']
score_files_dir  = pars['scoredir']
dcut             = pars['dcut']

dcuts = [dcut]
d_cuts.sort()

# read cliques
Data = pickle.load(open(clique_fname))
Residues, Cliques = [Data['Residues'], Data['Cliques']]
order = Cliques[0]['order']

# read different set of score files
Scores_lib = read_Scores_lib(score_files_dir, d_cuts, order)
depth_levels = []
for d_cut in d_cuts:
	# read residue clique score at different cut-off distances
	score_fname = score_files_dir + '/clique_score.'+str(order)+'-'+str(d_cut)+'.dat'
	scores_lib = read_table(score_fname, header = True)
	_depth_levels = scores_lib['DEPTH']
	scores_lib = dict(zip(zip(scores_lib['CLIQUE'], scores_lib['DEPTH']), scores_lib['LOG_ODD']))
	
	Scores_lib[d_cut] = scores_lib
	if depth_levels == []:
		depth_levels = set(_depth_levels)
	else:
		if depth_levels != set(_depth_levels):
			print 'error: different score files have different depth levels. Exiting'
			sys.exit(0)
		# end if
	# end if
# end for

depth_levels = Scores_lib['depth levels']
depth_cuts = []
for d in depth_levels:
	cut = [float(t) for t in d.split(',')]
	depth_cuts.append(cut)
# end for


# score cliques
for clique in Cliques:
	clique = score_clique(clique, d_cuts, Scores_lib)
# end for

# print output
# - use original filename is alright (only add information, no modification / deletion)
pickle.dump(Data, open(clique_fname, 'w'), pickle.HIGHEST_PROTOCOL)
