# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.


from msa_entropy import *
import sys

# input file are blocks of aligned characters
fname = sys.argv[1]
outfname = sys.argv[2]
M = open(fname).read().strip().split('\n')

# core function call
H, p = MSA_entropy(M) # H is JS-divergence, p is probability of each amino-acid type

# print output
fout = open(outfname, 'w')
aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
header = ['Pos', ' JSE'] + ['P('+i+')' for i in aa_list]
header_line = '\t'.join(header)+'\n'
fout.writelines(header_line)
for i in range(len(H)):
	output = [str(i+1), str(round(H[i], 4))] + [str(round(p[i][a], 4)) for a in range(len(aa_list))]
	outline = '\t'.join(output)+'\n'
	fout.writelines(outline)
# end for
fout.close()
