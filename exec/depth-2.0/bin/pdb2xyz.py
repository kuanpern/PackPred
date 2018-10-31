# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# from pdb coordinate to xyz coordinate
from PDB import *
import sys

def atom_name(atom):
	return ':'.join([atom.name, atom.residue.num, atom.residue.chain.name])
# end def

infile = sys.argv[1]
outfname = infile + '.xyz'
fout = open(outfname, 'w')
mdl = model(sys.argv[1], keywords = ['ATOM', 'HETATM'], remove_H = False, remove_alt = True, remove_non_std = False)
N = len(mdl.atoms)
for atom in mdl.atoms:
	out = '\t'.join([atom_name(atom), str(round(atom.x, 3)), str(round(atom.y, 3)), str(round(atom.z, 3))])
	fout.writelines(out+'\n')
# end for
fout.close()
