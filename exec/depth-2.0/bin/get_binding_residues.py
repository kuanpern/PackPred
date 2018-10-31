# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

import sys
from PDB import *

fname = sys.argv[1]
mdl = PDB(fname)

binding_residues = []
for i in range(len(mdl)):
	if mdl.T(i) == 1:
		res = str(mdl.resSeq(i))+':'+mdl.chainID(i)
		if res not in binding_residues:
			binding_residues.append(res)
		# end if
	# end if
# end for

out = ' '.join(binding_residues)
print out
