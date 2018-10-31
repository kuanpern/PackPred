# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# this script combines entropy and depth predicion values, and gives different pdb files for different combinations
from numpy import *
from PDB import *
import sys
pred_fname    = sys.argv[1]
entropy_fname = sys.argv[2]

def mean_float(s):
	v = []
	for t in s:
		try:
			v.append(float(t))
		except ValueError:
			continue
		# end try
	# end for
	return mean(v)
# end def 

def find_quartiles(S):
	mid = median(S)
	lower = [t for t in S if t < mid]
	upper = [t for t in S if t > mid]
	return median(lower), mid, median(upper)
# end def

root = pred_fname.split('/')[-1]

# 1. check if length of pdb file and entropy file match
entropy_data = open(entropy_fname).read().strip().split('\n')[1:]
entropy_data = [line.split() for line in entropy_data]
mdl = PDB(pred_fname)

# 1.1 read depth prediction
prediction = [mdl.T(i) for i in range(len(mdl)) if mdl.name(i) == 'CA']
# 1.2 read entropies <- 
# entropies = [float(t[1]) for t in entropy_data]
# following 10 lines make up for NA value in entropies (for benchmarking set only <- with unannotated residues)
buffer_entropy = [t[1] for t in entropy_data]
u = mean_float(buffer_entropy)
entropies = []
for t in buffer_entropy:
	try:
		entropies.append(float(t))
	except ValueError:
		entropies.append(u)
	# end try
# end for

# normalization to prevent 0 (error in logarithm function)
prediction = [max(t, 0) for t in prediction]
s_initial = sum(prediction)
prediction = array(prediction)+0.01
s_normalized = sum(prediction)
ratio = s_initial / s_normalized 
prediction = list(prediction * ratio)

# normalization of entropies (such that the distribution is similar to depth prediction)
#quartiles_entropy    = find_quartiles( entropies)
#quartiles_prediction = find_quartiles(prediction)
#median_entropy   , IQR_entropy    = [quartiles_entropy   [1], quartiles_entropy   [2] - quartiles_entropy   [0]]
#median_prediction, IQR_prediction = [quartiles_prediction[1], quartiles_prediction[2] - quartiles_prediction[0]]
#entropies = [x - median_entropy) / IQR_entropy * IQR_prediction + median_prediction for x in entropies]

adj_ratio = ( max(prediction) - min(prediction) ) / ( max(entropies) - min(entropies) )

entropies = array(entropies)
entropies = (entropies - min(entropies)) * adj_ratio + min(prediction)

#entropies = [min(x, 1) for x in entropies]
#entropies = [max(x, 0) for x in entropies]
#s_initial = sum(entropies)
#entropies = array(entropies)+0.01
#s_normalized = sum(entropies)
#ratio = s_initial / s_normalized 
entropies = list(entropies)

if len(entropies) != len(prediction):
	print 'data length unmatched', pred_fname, '(', len(prediction),')', entropy_fname, '(', len(entropies), ')'
	sys.exit(1)
# end if


def assign_T(mdl, p):
	residue_name_last = None
	c = -1
	for i in range(len(mdl)):
		residue_name = str(int(mdl.resSeq(i))) + ':' + mdl.chainID(i)
		if residue_name != residue_name_last:
			residue_name_last = residue_name
			c = c + 1
		# end if
		mdl._PDB__T[i] = p[c]
	# end for
	return mdl
# end def

# try simple combination here
c1 = linspace(0, 1, 21)
c2 = 1 - c1

for i in range(len(c1)):
	z = []
	for j in range(len(entropies)):
		v = 0.5*(c1[i]*float(prediction[j]) + c2[i]*float(entropies[j]))
		z.append(v)
	# end for
	mdl = assign_T(mdl, z)

	label = str(c1[i]) + '-' + str(c2[i])
	mdl.write(root+'.'+label+'.pdb')
# end for

