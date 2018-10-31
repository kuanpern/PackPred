# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# this script extract fasta sequence from a PDB file
from numpy import *
from PDB import *
import sys

code = {"A":"ALA", "R":"ARG", "D":"ASP", "N":"ASN", "C":"CYS", "E":"GLU", "Q":"GLN", "G":"GLY", "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS", "M":"MET", "F":"PHE", "P":"PRO", "S":"SER", "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL","ALA":"A", "ARG":"R", "ASP":"D", "ASN":"N", "CYS":"C", "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}

def extract_sequence(mdl):
	seq = {}
	residues = []
	for i in range(len(mdl)):
		c = mdl.chainID(i)
		residue_id = str(mdl.resSeq(i))+':'+c
		if residue_id not in residues:
			residues.append(residue_id)
			try:
				seq[c] = seq[c] + code[mdl.resName(i)]
			except KeyError:
				seq.update({c:''})
				seq[c] = seq[c] + code[mdl.resName(i)]
			# end try
		# end if
	# end for

	return seq
# end def


def fasta_format(seq, title, char_n = 70):
	out = '>'
	out = out + title + ';\n'

	# the sequence
	cuts = arange(0, len(seq), char_n)
	content = []
	for i in range(len(cuts)-1):
		content.append(seq[cuts[i]:cuts[i+1]])
	# end for
	content.append(seq[cuts[-1]:])
	content = '\n'.join(content).strip()+'\n'

	out = out + content
	return out
# end def

# read input
fname    = sys.argv[1]
out_root = sys.argv[2]

# get sequence
mdl = PDB(fname)
try:
	mdl.write("tmp.pdb")
except:
	pass
seq = extract_sequence(mdl)

# write output
chains = seq.keys()
for chain in chains:
	fasta_lines = fasta_format(seq[chain], out_root+'_'+chain)

	outfile = out_root+'_'+chain+'.fasta'
	print outfile
	fout = open(outfile, 'w')
	fout.writelines(fasta_lines+'\n')
	fout.close()
# end for
